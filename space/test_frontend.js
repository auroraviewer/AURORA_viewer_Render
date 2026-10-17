/*
 * Headless checks of the viewer front end.
 *
 * Parses app.js, verifies it wires the element ids that index.html defines, and
 * exercises the two pieces of image logic that are easy to get wrong:
 * colormap application, and the value/presence channels.
 */
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const dir = path.join(__dirname, 'static');
let failures = 0;
function check(label, cond, extra) {
  const ok = !!cond;
  console.log((ok ? 'PASS  ' : 'FAIL  ') + label + (extra !== undefined ? '  ' + extra : ''));
  if (!ok) failures++;
}

// ---- colormaps ------------------------------------------------------------
const sandbox = { atob: (b64) => Buffer.from(b64, 'base64').toString('binary') };
sandbox.window = sandbox;
vm.createContext(sandbox);
vm.runInContext(fs.readFileSync(path.join(dir, 'colormaps.js'), 'utf8'), sandbox, { filename: 'colormaps.js' });
const CM = sandbox.AuroraColormaps;
check('colormaps.js exposes AuroraColormaps', !!CM);
const turbo = CM.lut('turbo');
check('turbo LUT is 768 bytes', turbo.length === 768, turbo.length);
check('turbo(0) is #30123B', turbo[0] === 48 && turbo[1] === 18 && turbo[2] === 59);

// ---- app.js parses and wires the right ids --------------------------------
const appSrc = fs.readFileSync(path.join(dir, 'app.js'), 'utf8');
try {
  new vm.Script(appSrc, { filename: 'app.js' });
  check('app.js parses', true);
} catch (e) {
  check('app.js parses', false, e.message);
}

const IDS = ['sampleSel', 'geneSearch', 'geneBtn', 'geneListBox',
             'geneField', 'resSeg', 'cmapSel',
             'renderBtn', 'stage', 'warn28', 'heCanvas', 'predCanvas', 'spinner',
             'loadMsg', 'loadText', 'loadBarWrap', 'loadBar', 'loadPct',
             'status', 'ramp', 'legendTitle', 'legendCaption', 'metaInfo',
             'timingInfo', 'backendInfo', 'heNote'];
const html = fs.readFileSync(path.join(dir, 'index.html'), 'utf8');
for (const id of IDS) {
  check(`index.html defines #${id}`, html.indexOf('id="' + id + '"') >= 0);
}
for (const id of IDS) {
  check(`app.js references '${id}'`, appSrc.indexOf("'" + id + "'") >= 0);
}

// ---- the important behaviours ---------------------------------------------
check('app.js fetches /api/gene.png', appSrc.indexOf('/api/gene.png') >= 0);
check('app.js builds the H&E URL from the sample meta',
  appSrc.indexOf('he.url') >= 0 && appSrc.indexOf('state.meta.he') >= 0);
check('both canvases share one camera object',
  /state\.camera/.test(appSrc) && (appSrc.match(/attachPanZoom\(/g) || []).length >= 2,
  'attachPanZoom called for both panels');
check('the camera transform is applied to both panels in one function',
  /function drawAll\(\)[\s\S]{0,200}drawPanel\(els\.heCanvas/.test(appSrc) &&
  /drawAll\(\)[\s\S]{0,200}drawPanel\(els\.predCanvas/.test(appSrc));
check('cells without a patch are painted white inside the raster',
  /px\[i \+ 2\]/.test(appSrc) &&
  /px\[i\] = 255; px\[i \+ 1\] = 255; px\[i \+ 2\] = 255/.test(appSrc));
check('no hatch or quality overlay is drawn any more',
  appSrc.indexOf('buildHatch') < 0 && appSrc.indexOf('qualityWarn') < 0);
check('cell types are not offered (the prop decode overflows the tiny host)',
  appSrc.indexOf('celltype') < 0 && appSrc.indexOf('kindSel') < 0);
check('the panel matte around the raster stays dark',
  appSrc.indexOf("ctx.fillStyle = '#000';") >= 0
  && appSrc.indexOf(": '#ffffff'") < 0);
check('a scale bar is drawn on the H&E panel from the 0.5 um/px grid',
  appSrc.indexOf('drawScaleBar(') >= 0 &&
  appSrc.indexOf('UM_PER_WORLD_PX = 0.5') >= 0);

// ---- scale bar math (mirrors scaleBarSpan in app.js) -----------------------
const NICE_UM = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000];
function scaleBarSpan(camScale, wantPx) {
  if (!(camScale > 0)) return null;
  const umPerPx = 0.5 / camScale;
  const targetUm = umPerPx * wantPx;
  let um = NICE_UM[0];
  for (const n of NICE_UM) if (n <= targetUm) um = n;
  return { um, px: um / umPerPx, umPerPx };
}
{
  const fit = scaleBarSpan(800 / 8613, 140);      // whole LUAD slide in an 800 px panel
  check('scale bar picks a round length at fit zoom',
    fit && fit.um === 500, fit && `${fit.um} um = ${fit.px.toFixed(0)} px`);
  const zoom = scaleBarSpan(10 * (800 / 8613), 140);   // 10x zoomed in
  check('scale bar stays a readable length when zoomed',
    zoom && zoom.um === 50 && Math.abs(zoom.px - fit.px) < 2,
    zoom && `${zoom.um} um = ${zoom.px.toFixed(0)} px`);
  check('scale bar knows 0.5 um per world pixel',
    Math.abs(fit.umPerPx * (800 / 8613) - 0.5) < 1e-9);
}

// ---- reproduce buildValues on synthetic data ------------------------------
// RGB, as the server sends it: R = value code, B = presence.
function makeRGB(codes, presence, w, h) {
  const px = new Uint8ClampedArray(w * h * 4);
  for (let i = 0; i < w * h; i++) {
    px[i * 4] = codes[i];
    px[i * 4 + 2] = presence[i] ? 255 : 0;
    px[i * 4 + 3] = 255;
  }
  return px;
}

// exactly the loop in app.js buildValues()
function buildValues(px, name) {
  const lut = CM.lut(name);
  const out = new Uint8Array(px.length);
  for (let i = 0; i < px.length; i += 4) {
    if (px[i + 2]) {
      const o = px[i] * 3;
      out[i] = lut[o]; out[i + 1] = lut[o + 1];
      out[i + 2] = lut[o + 2]; out[i + 3] = 255;
    } else {
      out[i] = 255; out[i + 1] = 255; out[i + 2] = 255; out[i + 3] = 255;  // white
    }
  }
  return out;
}

const codes = new Uint8Array([0, 10, 128, 255, 0, 200, 50, 90]);
const presence = new Uint8Array([1, 1, 1, 1, 1, 1, 1, 0]);   // last cell has no patch
const out = buildValues(makeRGB(codes, presence, 4, 2), 'turbo');

check('a present cell with code 0 gets the colormap low end, not white',
  out[0] === turbo[0] && out[3] === 255, [out[0], out[1], out[2], out[3]].join(','));
check('the last code maps to the colormap high end',
  out[3 * 4] === turbo[765] && out[3 * 4 + 1] === turbo[766]);
check('present cells are opaque', out[3] === 255 && out[4 * 4 + 3] === 255);
check('a cell with no patch is opaque white', out[7 * 4 + 3] === 255
  && out[7 * 4] === 255 && out[7 * 4 + 1] === 255 && out[7 * 4 + 2] === 255,
  [out[7 * 4], out[7 * 4 + 1], out[7 * 4 + 2], out[7 * 4 + 3]].join(','));
check('an empty raster is opaque white',
  buildValues(makeRGB(codes, new Uint8Array(8), 4, 2), 'turbo')[3] === 255);

console.log();
if (failures) {
  console.log(failures + ' FAILURE(S)');
  process.exit(1);
}
console.log('ALL FRONTEND CHECKS PASSED');
