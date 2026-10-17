/**
 * Headless test of the sample-loading UI.
 *
 * Pins three behaviours the user asked for:
 *   1. the warning message stays visible while a sample downloads
 *   2. the progress bar reflects the server's per-file progress
 *   3. switching to another sample shows the bar again
 *
 * The server's prefetch endpoint is scripted, so the test is deterministic.
 */
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const DIR = path.join(__dirname, 'static');
const appSrc = fs.readFileSync(path.join(DIR, 'app.js'), 'utf8');
const cmSrc = fs.readFileSync(path.join(DIR, 'colormaps.js'), 'utf8');
const html = fs.readFileSync(path.join(DIR, 'index.html'), 'utf8');

const failures = [];
function check(label, ok, extra) {
  console.log((ok ? 'PASS  ' : 'FAIL  ') + label + (extra !== undefined ? `  ${extra}` : ''));
  if (!ok) failures.push(label);
}

// ---- mock DOM (same shape as test_boot.js) --------------------------------
const listeners = new Map();
function makeEl(id) {
  const el = {
    id, style: {}, dataset: {}, className: '', textContent: '', innerHTML: '',
    value: '', width: 800, height: 600, files: [], children: [],
    classList: { add() {}, remove() {}, toggle() {}, contains: () => false },
    setAttribute() {}, getAttribute: () => null, removeAttribute() {},
    appendChild(c) { el.children.push(c); return c; },
    removeChild() {}, focus() {}, blur() {}, click() {},
    setPointerCapture() {}, releasePointerCapture() {},
    addEventListener(type, fn) {
      const k = `${id}:${type}`;
      if (!listeners.has(k)) listeners.set(k, []);
      listeners.get(k).push(fn);
    },
    removeEventListener() {},
    getBoundingClientRect: () => ({ left: 0, top: 0, width: 800, height: 600 }),
    getContext() {
      return {
        canvas: el,
        setTransform() {}, resetTransform() {}, fillRect() {}, clearRect() {},
        drawImage() {}, fillText() {}, beginPath() {}, closePath() {}, moveTo() {},
        lineTo() {}, arc() {}, stroke() {}, fill() {}, save() {}, restore() {},
        clip() {}, rect() {}, translate() {}, scale() {}, rotate() {}, strokeRect() {},
        getImageData: (x, y, w, h) => ({ data: new Uint8ClampedArray(Math.max(4, w * h * 4)) }),
        createImageData: (w, h) => ({ data: new Uint8ClampedArray(Math.max(4, w * h * 4)) }),
        putImageData() {}, createPattern: () => ({}),
        measureText: (t) => ({ width: (t || '').length * 6 }),
        set imageSmoothingEnabled(v) {}, get imageSmoothingEnabled() { return true; }
      };
    }
  };
  return el;
}
const defined = new Set([...html.matchAll(/id="([^"]+)"/g)].map((m) => m[1]));
const els = {};
for (const id of defined) els[id] = makeEl(id);

const documentMock = {
  readyState: 'complete',
  activeElement: null,
  documentElement: makeEl('html'),
  getElementById: (id) => els[id] || null,
  createElement: (t) => makeEl(t),
  addEventListener(type, fn) {
    const k = `document:${type}`;
    if (!listeners.has(k)) listeners.set(k, []);
    listeners.get(k).push(fn);
  },
  querySelector: () => null,
  querySelectorAll: () => [],
  body: makeEl('body')
};

// ---- scripted server ------------------------------------------------------
// mode: "download" -> prefetch reports progress; "done" -> finished
const job = { mode: 'download', calls: 0 };
const META = {
  genes: ['EPCAM', 'TUMOR'],
  cell_types: [],
  layers: [{ cell_px: 14, microns: 7, grid: [4, 3], extent: [40, 30],
             origin: [0, 0], patches: 12 }],
  he: { available: false }
};
function jsonFor(url) {
  if (url.indexOf('/prefetch') >= 0) {
    job.calls++;
    if (job.mode === 'done') return { state: 'done', file_progress: 1 };
    const second = job.calls > 1;         // pretend it keeps going
    return {
      state: 'running',
      current: 'SAMPLE-highres-pred.h5',
      files_total: 2, files_done: 0,
      file_downloaded: second ? 800e6 : 420e6,
      file_total: 1000e6,
      file_progress: second ? 0.8 : 0.42
    };
  }
  if (url.indexOf('/meta') >= 0) return META;
  if (url.indexOf('/api/samples') >= 0) {
    return {
      count: 2,
      samples: [{ sample_id: 'S1', size: 10 }, { sample_id: 'S2', size: 20 }],
      default: 'S1', default_gene: 'EPCAM'
    };
  }
  if (url.indexOf('/api/health') >= 0) return { backend: 'hf', he_images: { indexed: 0 } };
  return {};
}

const sandbox = {
  console,
  document: documentMock,
  fetch: (url) => Promise.resolve({
    ok: true, status: 200,
    json: () => Promise.resolve(jsonFor(url)),
    text: () => Promise.resolve(JSON.stringify(jsonFor(url)))
  }),
  setTimeout, clearTimeout, setInterval, clearInterval,
  requestAnimationFrame: (fn) => setTimeout(fn, 0),
  atob: (s) => Buffer.from(s, 'base64').toString('binary'),
  btoa: (s) => Buffer.from(s, 'binary').toString('base64'),
  performance: { now: () => 0 },
  Image: function () { return makeEl('img'); },
  Path2D: function () { return { addPath() {}, moveTo() {}, lineTo() {}, closePath() {} }; },
  location: { hash: '', href: 'http://127.0.0.1:7860/' },
  navigator: { userAgent: 'node' },
  URLSearchParams,
  parseInt, parseFloat, isNaN, Math, Date, JSON
};
sandbox.addEventListener = function () {};
sandbox.removeEventListener = function () {};
sandbox.window = sandbox;
sandbox.globalThis = sandbox;
sandbox.self = sandbox;

const ctx = vm.createContext(sandbox);
vm.runInContext(cmSrc, ctx, { filename: 'colormaps.js' });
vm.runInContext(appSrc, ctx, { filename: 'app.js' });

const wait = (ms) => new Promise((r) => setTimeout(r, ms));

function barWidth() { return els.loadBar.style.width; }
function barShown() { return els.loadBarWrap.style.display === 'block'; }

(async () => {
  await wait(120);
  check('a download shows the progress bar', barShown(), barShown() ? 'shown' : 'hidden');
  check('the bar reflects the server progress', barWidth() === '42%', barWidth());
  const label = els.loadPct.textContent || '';
  check('the bar names the file and percentage',
    label.indexOf('SAMPLE-highres-pred.h5') >= 0 && label.indexOf('42%') >= 0, label);
  check('the warning message is still visible during the download',
    els.loadMsg.style.display === 'block' &&
    (els.loadText.textContent || '').indexOf('Loading sample data') >= 0,
    els.loadText.textContent);

  await wait(900);                       // a second poll
  check('the bar keeps moving as bytes arrive', barWidth() === '80%', barWidth());

  job.mode = 'done';                     // finish the download
  await wait(900);
  check('the bar is hidden once the sample is local', !barShown(),
    els.loadBarWrap.style.display || '(unset)');

  // ---- switch to another sample ------------------------------------------
  job.mode = 'download';
  job.calls = 0;
  els.sampleSel.value = 'S2';
  const onChange = (listeners.get('sampleSel:change') || [])[0];
  check('the sample picker is wired', typeof onChange === 'function');
  onChange();

  await wait(120);
  check('switching samples shows the progress bar again', barShown(),
    barShown() ? 'shown' : 'hidden');
  check('the bar restarts for the new sample', barWidth() === '42%', barWidth());
  check('the warning message is visible on a switch too',
    els.loadMsg.style.display === 'block' &&
    (els.loadText.textContent || '').indexOf('Loading sample data') >= 0,
    els.loadText.textContent);

  job.mode = 'done';
  await wait(900);
  check('the bar hides again after the switch finishes', !barShown());

  console.log();
  if (failures.length) {
    console.log(`${failures.length} FAILURE(S): ` + failures.join('; '));
    process.exit(1);
  }
  console.log('LOADING UI CHECKS PASSED');
  process.exit(0);
})();
