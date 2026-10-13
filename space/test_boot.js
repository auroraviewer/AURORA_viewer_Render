/**
 * Execute the real static/app.js against a mock DOM and assert that start-up
 * actually succeeds.
 *
 * This exists because a broken app.js fails silently: the page still loads, the
 * server logs look fine, and the only symptom is that no /api/samples request is
 * ever made.  Running the file is the only way to catch that class of bug.
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

// ---- ids ------------------------------------------------------------------
const defined = new Set([...html.matchAll(/id="([^"]+)"/g)].map((m) => m[1]));
const want = [...appSrc.matchAll(/var ids = \[([\s\S]*?)\];/g)]
  .flatMap((m) => [...m[1].matchAll(/'([^']+)'/g)].map((x) => x[1]));
const missingIds = want.filter((id) => !defined.has(id));
check('every id app.js asks for exists in index.html',
  missingIds.length === 0, missingIds.join(', ') || 'all present');

// ---- mock DOM -------------------------------------------------------------
const listeners = new Map();
const calls = [];

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
        clip() {}, rect() {}, translate() {}, scale() {}, rotate() {},
        getImageData: (x, y, w, h) => ({ data: new Uint8ClampedArray(Math.max(4, w * h * 4)) }),
        createImageData: (w, h) => ({ data: new Uint8ClampedArray(Math.max(4, w * h * 4)) }),
        putImageData() {}, createPattern: () => ({}),
        measureText: () => ({ width: 0 }),
        set imageSmoothingEnabled(v) {}, get imageSmoothingEnabled() { return true; }
      };
    }
  };
  return el;
}

const store = {};
for (const id of defined) store[id] = makeEl(id);

const documentMock = {
  readyState: 'complete',
  activeElement: null,
  documentElement: makeEl('html'),
  getElementById: (id) => store[id] || null,
  createElement: (t) => makeEl(t),
  addEventListener(type, fn) { (listeners.get(`document:${type}`) || listeners.set(`document:${type}`, []).get(`document:${type}`)).push(fn); },
  querySelector: () => null,
  querySelectorAll: () => [],
  body: makeEl('body')
};

const jsonFor = (url) => {
  if (url.includes('/api/samples')) {
    return { count: 1, samples: [{ sample_id: 'S1', size: 10 }] };
  }
  if (url.includes('/api/health')) {
    return { backend: 'local', he_images: { indexed: 2 }, samples_loaded: 1 };
  }
  return { count: 0, samples: [] };
};

const sandbox = {
  console,
  document: documentMock,
  fetch: (url) => {
    calls.push(url);
    return Promise.resolve({
      ok: true, status: 200,
      json: () => Promise.resolve(jsonFor(url)),
      text: () => Promise.resolve(JSON.stringify(jsonFor(url)))
    });
  },
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
sandbox.addEventListener = function (type, fn) {
  const k = 'window:' + type;
  if (!listeners.has(k)) listeners.set(k, []);
  listeners.get(k).push(fn);
};
sandbox.removeEventListener = function () {};
sandbox.window = sandbox;
sandbox.globalThis = sandbox;
sandbox.self = sandbox;

const ctx = vm.createContext(sandbox);

// ---- run it ---------------------------------------------------------------
let failure = null;
try {
  vm.runInContext(cmSrc, ctx, { filename: 'colormaps.js' });
  vm.runInContext(appSrc, ctx, { filename: 'app.js' });
} catch (e) {
  failure = e;
}

check('app.js runs without throwing', failure === null,
  failure ? `${failure.name}: ${failure.message}` : '');
if (failure) {
  console.log('      ' + String(failure.stack).split('\n')[1]);
}

setTimeout(() => {
  console.log();
  check('it requests /api/samples on start-up',
    calls.some((u) => u.includes('/api/samples')),
    calls.length ? calls.join(', ') : 'no requests at all');
  check('the page rendered without an unhandled error',
    listeners.size > 0, `${listeners.size} listeners attached`);

  console.log();
  if (failures.length) {
    console.log(`${failures.length} FAILURE(S): ` + failures.join('; '));
    process.exit(1);
  }
  console.log('FRONT END BOOTS CLEANLY');
}, 50);

