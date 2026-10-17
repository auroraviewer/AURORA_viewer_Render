/*
 * AURORA viewer front end.
 *
 * Two panels - H&E and prediction - share ONE camera (world scale + offset), so
 * panning or zooming either one moves both.  The prediction panel is rendered
 * from the server's grayscale PNG through a colormap LUT in the browser, which
 * is why changing gene or colormap never needs a new server-side artefact.
 *
 * The server sends an RGB PNG: R is the value code and B presence.  Presence 0
 * becomes transparent, so areas with no patch show the page through rather than
 * a rectangle of near-black.
 */
(function () {
  'use strict';

  var CM = window.AuroraColormaps;

  var state = {
    sample: null,
    gene: null,
    catalog: [],          // gene list loaded from LUAD_gene.csv (fallback: H5 genes)
    cell: 14,
    cmap: 'turbo',
    meta: null,
    layer: null,           // current layer descriptor
    bitmap: null,          // decoded prediction PNG
    values: null,          // offscreen canvas: colormapped RGB
    heBitmap: null,
    heInfo: null,
    world: { w: 1, h: 1 },
    camera: { scale: 1, tx: 0, ty: 0, fit: 1 }
  };

  var els = {};
  // Every id here must exist in index.html.  test_boot.js runs this file against
  // a mock DOM and fails if one is missing, or if start-up throws at all - the
  // failure mode being guarded against is a page that loads and then silently
  // does nothing.
  var ids = ['sampleSel', 'geneSearch', 'geneBtn', 'geneListBox',
             'geneField', 'resSeg', 'cmapSel', 'renderBtn', 'stage', 'warn28',
             'heCanvas', 'predCanvas', 'spinner', 'loadMsg', 'status',
             'ramp', 'legendTitle', 'legendCaption', 'metaInfo', 'timingInfo',
             'backendInfo', 'heNote'];

  function $(id) { return document.getElementById(id); }
  function api(p) {
    return fetch(p).then(function (r) {
      if (!r.ok) {
        return r.text().then(function (t) {
          var m = t; try { m = JSON.parse(t).detail || t; } catch (e) {}
          throw new Error(m || ('HTTP ' + r.status));
        });
      }
      return r.json();
    });
  }
  function busy(on, msg) {
    els.spinner.style.display = on ? 'block' : 'none';
    if (on && msg) {
      els.loadMsg.textContent = msg;
      els.loadMsg.style.display = 'block';
    } else {
      els.loadMsg.style.display = 'none';
    }
  }
  // ---- gene catalog -----------------------------------------------------
  function parseGeneCSV(text) {
    var out = [];
    (text || '').replace(/^\uFEFF/, '').split(/\r?\n/).forEach(function (line) {
      var parts = line.split(',');
      if (parts.length < 2) return;
      var g = parts[1].trim();
      if (!g || g === 'Genes') return;
      out.push(g);
    });
    return out;
  }
  function catalogGenes() {
    if (state.catalog && state.catalog.length) return state.catalog;
    return state.meta ? state.meta.genes : [];
  }
  function inSample(g) {
    return !!(state.meta && state.meta.genes.indexOf(g) >= 0);
  }
  function availGenes() {
    var out = [];
    catalogGenes().forEach(function (g) { if (inSample(g)) out.push(g); });
    if (out.length) return out;
    return state.meta ? state.meta.genes : [];
  }
  function status(msg, err) {
    if (!msg) { els.status.style.display = 'none'; return; }
    els.status.textContent = msg;
    els.status.className = err ? 'error' : '';
    els.status.style.display = 'block';
  }

  // ------------------------------------------------------------ camera
  function panelSize(cv) {
    var r = cv.parentElement.getBoundingClientRect();
    var dpr = window.devicePixelRatio || 1;
    var w = Math.max(1, Math.round(r.width)), h = Math.max(1, Math.round(r.height));
    if (cv.width !== w * dpr || cv.height !== h * dpr) {
      cv.width = w * dpr; cv.height = h * dpr;
      cv.style.width = r.width + 'px'; cv.style.height = r.height + 'px';
    }
    return { w: w, h: h, dpr: dpr };
  }

  function fitCamera() {
    var cam = state.camera;
    var s = panelSize(els.predCanvas);
    var sx = s.w / state.world.w, sy = s.h / state.world.h;
    cam.fit = Math.min(sx, sy) * 0.94;
    cam.scale = cam.fit;
    cam.tx = (s.w - state.world.w * cam.scale) / 2;
    cam.ty = (s.h - state.world.h * cam.scale) / 2;
  }

  function clampCamera() {
    var cam = state.camera;
    var s = panelSize(els.predCanvas);
    var m = 60;
    cam.tx = Math.min(s.w - m, Math.max(m - state.world.w * cam.scale, cam.tx));
    cam.ty = Math.min(s.h - m, Math.max(m - state.world.h * cam.scale, cam.ty));
  }

  function drawAll() {
    drawPanel(els.heCanvas, state.heBitmap, true);
    drawPanel(els.predCanvas, state.values, false);
  }

  function drawPanel(cv, bmp, isHE) {
    var s = panelSize(cv);
    var ctx = cv.getContext('2d');
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.imageSmoothingEnabled = state.camera.scale < 1;
    ctx.fillStyle = '#000';
    ctx.fillRect(0, 0, cv.width, cv.height);

    if (isHE && !bmp) {
      ctx.setTransform(s.dpr, 0, 0, s.dpr, 0, 0);
      ctx.fillStyle = '#1b1e26';
      ctx.fillRect(0, 0, s.w, s.h);
      ctx.fillStyle = '#7d8798';
      ctx.font = '13px system-ui, sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText('No H&E image for this sample', s.w / 2, s.h / 2 - 8);
      ctx.fillText('set AURORA_HE_DIR to enable the side-by-side view',
                   s.w / 2, s.h / 2 + 12);
      return;
    }
    if (!bmp) return;

    var cam = state.camera;
    ctx.setTransform(s.dpr, 0, 0, s.dpr, 0, 0);
    ctx.globalAlpha = 1;

    if (isHE) {
      // The H&E sits at the world origin and spans its native pixel size.
      // It is NOT stretched over the whole world: the prediction may be narrower
      // or shift into negative coordinates, and stretching would misalign them.
      var hw = state.heInfo.width * cam.scale;
      var hh = state.heInfo.height * cam.scale;
      ctx.drawImage(bmp, cam.tx, cam.ty, hw, hh);
    } else if (state.predRect) {
      // The prediction is placed at its reported origin and extent, so a layer
      // starting at a negative coordinate or narrower than the image lands
      // exactly where it belongs.
      var r = state.predRect;
      var px = cam.tx + r.x * cam.scale;
      var py = cam.ty + r.y * cam.scale;
      var pw = r.w * cam.scale;
      var ph = r.h * cam.scale;
      // Very coarse levels (112 um is a 197x45 raster) must NOT be smoothed when
      // scaled up, or neighbouring cells blend into diagonal smears.  Keep the
      // cells crisp so the resolution is honest.
      var coarse = bmp.width > 0 && (pw / bmp.width) >= 2.5;
      var smooth = cam.scale < 1 && !coarse;
      ctx.imageSmoothingEnabled = isHE ? (cam.scale < 1) : smooth;
      ctx.drawImage(bmp, px, py, pw, ph);
    }
    ctx.globalAlpha = 1;
  }

  // ---------------------------------------------------- interaction
  /**
   * Pan and zoom one panel by driving the SHARED camera, so the other panel
   * follows.  There is no per-panel transform: whatever is done here shows up in
   * both, which is the whole point of the side-by-side view.
   */
  function attachPanZoom(cv) {
    if (!cv) return;
    var host = cv.parentElement || cv;
    var dragging = false, lx = 0, ly = 0;

    host.addEventListener('pointerdown', function (e) {
      if (e.button !== undefined && e.button !== 0) return;
      dragging = true; lx = e.clientX; ly = e.clientY;
      if (host.setPointerCapture) { try { host.setPointerCapture(e.pointerId); } catch (err) {} }
      host.style.cursor = 'grabbing';
      e.preventDefault();
    });

    host.addEventListener('pointermove', function (e) {
      if (!dragging) return;
      state.camera.tx += e.clientX - lx;
      state.camera.ty += e.clientY - ly;
      lx = e.clientX; ly = e.clientY;
      clampCamera(); drawAll();
    });

    function endDrag() {
      if (!dragging) return;
      dragging = false;
      host.style.cursor = 'grab';
    }
    host.addEventListener('pointerup', endDrag);
    host.addEventListener('pointercancel', endDrag);
    host.addEventListener('pointerleave', endDrag);

    // zoom about the pointer, so the slide point under the cursor stays put
    host.addEventListener('wheel', function (e) {
      e.preventDefault();
      var s = panelSize(cv);
      var r = host.getBoundingClientRect ? host.getBoundingClientRect()
                                         : { left: 0, top: 0 };
      var mx = (e.clientX || 0) - r.left;
      var my = (e.clientY || 0) - r.top;
      var cam = state.camera;
      var f = Math.exp(-(e.deltaY || 0) * 0.0015);
      var next = Math.min(40, Math.max(0.02, cam.scale * f));
      var k = next / cam.scale;
      cam.tx = mx - (mx - cam.tx) * k;
      cam.ty = my - (my - cam.ty) * k;
      cam.scale = next;
      clampCamera(); drawAll();
    }, { passive: false });

    // double click fits the view again
    host.addEventListener('dblclick', function () {
      fitCamera(); drawAll();
    });
  }

  // ---------------------------------------------------- colour mapping
  function buildValues(bmp) {
    var c = document.createElement('canvas');
    c.width = bmp.width; c.height = bmp.height;
    var ctx = c.getContext('2d', { willReadFrequently: true });
    ctx.drawImage(bmp, 0, 0);
    var data = ctx.getImageData(0, 0, c.width, c.height);
    var px = data.data;
    var lut = CM.lut(state.cmap);

    // The server sends RGB: R is the value code, B is presence.
    //   * presence 0  -> alpha 0, so "no patch" is genuinely transparent and the
    //     panel shows through instead of a rectangle of near-black.
    //   * presence 255 -> the value is mapped through the colormap, so a real
    //     expression of 0 gets the colormap's low colour, not transparency.
    for (var i = 0; i < px.length; i += 4) {
      var present = px[i + 2];
      if (present) {
        var o = px[i] * 3;
        px[i] = lut[o]; px[i + 1] = lut[o + 1]; px[i + 2] = lut[o + 2];
        px[i + 3] = 255;
      } else {
        px[i + 3] = 0;
      }
    }
    ctx.putImageData(data, 0, 0);
    state.values = c;
  }


  // ---------------------------------------------------- loading
  function loadSamples() {
    return api('/api/samples').then(function (d) {
      els.sampleSel.innerHTML = '';
      if (!d.samples.length) {
        status('No samples found. Set AURORA_LOCAL_DIR or AURORA_HF_REPO.', true);
        return;
      }
      d.samples.forEach(function (s) {
        var o = document.createElement('option');
        o.value = s.sample_id; o.textContent = s.sample_id;
        els.sampleSel.appendChild(o);
      });
      var h = readHash();
      var want = (h.sample && d.samples.some(function (s) { return s.sample_id === h.sample; }))
        ? h.sample
        : (d.default && d.samples.some(function (s) { return s.sample_id === d.default; })
            ? d.default : d.samples[0].sample_id);
      els.sampleSel.value = want;
      return selectSample(want, h);
    }).catch(function (e) { status('Cannot list samples: ' + e.message, true); });
  }

  function selectSample(id, restore) {
    state.sample = id;
    busy(true, 'Loading sample data \u2014 the first load downloads it from '
      + 'Hugging Face and can take a minute or two. After that, switching genes '
      + 'and resolutions is fast.');
    status('Reading ' + id + '\u2026');
    return api('/api/samples/' + encodeURIComponent(id) + '/meta').then(function (m) {
      state.meta = m;
      els.geneListBox.style.display = 'none';

      // resolutions
      els.resSeg.innerHTML = '';
      m.layers.forEach(function (L) {
        var b = document.createElement('button');
        b.textContent = L.microns + ' \u00b5m';
        b.title = L.grid[0] + ' \u00d7 ' + L.grid[1] + ' cells'
        b.setAttribute('aria-pressed', String(L.cell_px === state.cell));
        b.onclick = function () { setCell(L.cell_px); };
        els.resSeg.appendChild(b);
      });
      if (!m.layers.some(function (L) { return L.cell_px === state.cell; })) {
        state.cell = m.layers[0].cell_px;
      }
      syncResButtons();
      updateWarn28();

      var gene = (restore && restore.gene && inSample(restore.gene))
        ? restore.gene : (availGenes()[0] || m.genes[0]);
      els.geneSearch.value = gene;
      state.gene = gene;
      if (restore && restore.cmap) { state.cmap = restore.cmap; els.cmapSel.value = restore.cmap; }

      setWorld();
      loadHE();
      busy(false); status('');
      return render(false);
    }).catch(function (e) {
      busy(false); status('Cannot read ' + id + ': ' + e.message, true);
    });
  }

  function setWorld() {
    // The world IS native image-pixel space: prediction coordinates are in
    // pixels of the rescaled image, and the H&E's native size is reported by the
    // server (the served image is downscaled, so its own width cannot be used).
    var L = currentLayer();
    var heW = state.heInfo ? state.heInfo.width : 0;
    var heH = state.heInfo ? state.heInfo.height : 0;

    var w = heW || (L ? L.extent[0] : 1);
    var h = heH || (L ? L.extent[1] : 1);
    // a layer may start at a negative coordinate or extend past the image
    if (L) {
      var o = L.origin || [0, 0];
      w = Math.max(w, o[0] + L.extent[0]);
      h = Math.max(h, o[1] + L.extent[1]);
    }
    state.world = { w: Math.max(1, w), h: Math.max(1, h) };
    updatePredRect();
    fitCamera();
  }

  // Where the current layer's raster sits, in native image pixels.
  function updatePredRect() {
    var L = currentLayer();
    if (!L) { state.predRect = null; return; }
    var o = (L && L.origin) || [0, 0];
    var e = (L && L.extent) || [1, 1];
    var ri = state.renderInfo;
    // Use the server's reported origin and extent.  The extent is measured on the
    // raster's coordinate axis, so it includes any hole in the patch lattice -
    // which native_grid * cell_px would silently close up, sliding everything
    // right of the hole across and butting two pieces of tissue together.
    if (ri && ri.origin && ri.extent) {
      state.predRect = { x: ri.origin[0], y: ri.origin[1], w: ri.extent[0], h: ri.extent[1] };
    } else {
      state.predRect = { x: o[0], y: o[1], w: e[0], h: e[1] };
    }
  }

  function currentLayer() {
    if (!state.meta) return null;
    for (var i = 0; i < state.meta.layers.length; i++) {
      if (state.meta.layers[i].cell_px === state.cell) return state.meta.layers[i];
    }
    return state.meta.layers[0] || null;
  }

  function loadHE() {
    var he = state.meta && state.meta.he;
    els.heNote.textContent = '';
    if (!he || !he.available) {
      state.heBitmap = null; state.heInfo = null;
      els.heNote.textContent = 'no H&E';
      drawAll();
      return;
    }
    if (he.placeholder) els.heNote.textContent = 'placeholder image';
    var img = new Image();
    img.onload = function () {
      state.heBitmap = img;
      // the served image is downscaled; the world is in native image pixels, so
      // keep the native size for placement and let the canvas scale the bitmap
      var native = he.image_size || [img.width, img.height];
      state.heInfo = { width: native[0], height: native[1] };
      setWorld();
      drawAll();
    };
    img.onerror = function () { state.heBitmap = null; state.heInfo = null; drawAll(); };
    img.src = he.url + '?max_width=6000';
  }

  function setCell(cell) {
    if (state.cell === cell) return;
    state.cell = cell;
    syncResButtons();
    updateWarn28();
    setWorld();
    render(true);
  }

  function updateWarn28() {
    els.warn28.style.display = (state.cell === 56) ? 'block' : 'none';
  }

  function syncResButtons() {
    Array.prototype.forEach.call(els.resSeg.children, function (b, i) {
      var L = state.meta.layers[i];
      b.setAttribute('aria-pressed', String(L && L.cell_px === state.cell));
    });
  }

  function render(keepView) {
    if (!state.sample) return;
    var name = state.gene;
    if (!name) return;
    var url = '/api/gene.png?sample=' + encodeURIComponent(state.sample)
        + '&gene=' + encodeURIComponent(name)
        + '&cell=' + state.cell + '&width=2080';

    busy(true);
    status('Rendering ' + name + ' at ' + (currentLayer() ? currentLayer().microns : '?') + ' \u00b5m\u2026');
    var t0 = performance.now();

    // ask where this view sits in image-pixel space, in parallel with the image
    var statsUrl = '/api/view/stats?sample=' + encodeURIComponent(state.sample)
      + '&name=' + encodeURIComponent(name)
      + '&kind=gene&cell=' + state.cell + '&width=2080';

    var img = new Image();
    var pending = 2;
    function done() {
      if (--pending > 0) return;
      if (!keepView) { setWorld(); } else { clampCamera(); }
      drawAll();
      busy(false); status('');
      els.timingInfo.textContent = 'rendered in ' + Math.round(performance.now() - t0) + ' ms';
      updateLegend();
      updateMetaLine();
      updateUrl();
    }

    api(statsUrl).then(function (st) {
      state.renderInfo = st;
      updatePredRect();
      done();
    }).catch(function (e) {
      state.renderInfo = null;
      status('Could not read view info: ' + e.message, true);
      done();
    });

    img.onload = function () {
      state.bitmap = img;
      buildValues(img);
      done();
    };
    img.onerror = function () {
      busy(false);
      status('Could not render ' + name + '.', true);
    };
    img.src = url;
  }

  function updateLegend() {
    CM.drawRamp(els.ramp, state.cmap);
    var name = state.gene;
    var ri = state.renderInfo;
    els.legendTitle.textContent = 'Predicted expression';
    var L = currentLayer();
    els.legendCaption.textContent = name
      + (L ? ' \u00b7 ' + L.microns + ' \u00b5m' : '');
  }

  function updateMetaLine() {
    var L = currentLayer();
    if (!L || !state.meta) { els.metaInfo.textContent = '\u00a0'; return; }
    els.metaInfo.textContent = state.sample + ' \u00b7 '
      + L.grid[0] + '\u00d7' + L.grid[1] + ' cells \u00b7 '
      + L.patches.toLocaleString() + ' patches \u00b7 '
      + L.microns + ' \u00b5m \u00b7 ' + L.extent[0] + '\u00d7' + L.extent[1] + ' px world';
  }

  // ---------------------------------------------------- url
  function updateUrl() {
    var h = '#sample=' + encodeURIComponent(state.sample)
      + '&gene=' + encodeURIComponent(state.gene)
      + '&cell=' + state.cell + '&cmap=' + state.cmap;
    if (history.replaceState) history.replaceState(null, '', h);
  }
  function readHash() {
    var out = {};
    (location.hash || '').replace(/^#/, '').split('&').forEach(function (kv) {
      if (!kv) return;
      var i = kv.indexOf('=');
      if (i < 0) return;
      out[decodeURIComponent(kv.slice(0, i))] = decodeURIComponent(kv.slice(i + 1));
    });
    if (out.cell) state.cell = parseInt(out.cell, 10) || 14;
    if (out.cmap) state.cmap = out.cmap;
    return out;
  }

  // ---------------------------------------------------- wiring
  function wire() {
    els.sampleSel.addEventListener('change', function () { selectSample(els.sampleSel.value, null); });
    els.cmapSel.addEventListener('change', function () {
      state.cmap = els.cmapSel.value;
      if (state.bitmap) { buildValues(state.bitmap); drawAll(); }
      updateLegend(); updateUrl();
    });
    els.renderBtn.addEventListener('click', function () { setWorld(); render(false); });

    // ---- gene picker: type to search, or open the full list ----
    function geneListVisible() { return els.geneListBox.style.display === 'block'; }
    function hideGeneList() { els.geneListBox.style.display = 'none'; }
    function pickGene(g) {
      if (!inSample(g)) {
        status('Gene "' + g + '" is not present in this sample\u2019s predictions.', true);
        return;
      }
      state.gene = g;
      els.geneSearch.value = g;
      hideGeneList();
      status('');
      render(true);
    }
    function renderGeneList(filter) {
      var genes = catalogGenes();
      var q = (filter || '').trim().toLowerCase();
      var filtering = q.length > 0;
      var list = filtering ? genes.filter(function (g) { return g.toLowerCase().indexOf(q) >= 0; })
                           : genes;
      var MAX = 300;
      els.geneListBox.innerHTML = '';
      var frag = document.createDocumentFragment();
      if (filtering) {
        // a filtered list may collapse to one gene; always offer the full list
        var all = document.createElement('div');
        all.className = 'geneItem all';
        all.textContent = 'Show all genes';
        all.addEventListener('click', function (e) {
          e.stopPropagation();          // the rebuild detaches this node; keep
          renderGeneList('');           // the list open after it
        });
        frag.appendChild(all);
      }
      for (var i = 0; i < Math.min(list.length, MAX); i++) {
        (function (g) {
          var ok = inSample(g);
          var d = document.createElement('div');
          d.className = ok ? 'geneItem' : 'geneItem off';
          d.title = ok ? '' : 'Not present in this sample\u2019s predictions';
          d.textContent = g;
          d.addEventListener('click', function () { pickGene(g); });
          frag.appendChild(d);
        })(list[i]);
      }
      if (list.length > MAX) {
        var more = document.createElement('div');
        more.className = 'geneItem';
        more.textContent = '\u2026 ' + (list.length - MAX) + ' more \u2014 keep typing to narrow';
        frag.appendChild(more);
      }
      els.geneListBox.appendChild(frag);
    }
    els.geneBtn.addEventListener('click', function () {
      if (geneListVisible()) { hideGeneList(); return; }
      renderGeneList(els.geneSearch.value);
      els.geneListBox.style.display = 'block';
    });
    document.addEventListener('click', function (e) {
      if (geneListVisible()
          && e.target !== els.geneBtn
          && e.target !== els.geneSearch
          && !els.geneListBox.contains(e.target)) {
        hideGeneList();
      }
    });

    var timer = null;
    function commitGene() {
      var v = els.geneSearch.value.trim();
      if (!v || !state.meta) return;
      if (!inSample(v)) {
        status('Gene "' + v + '" is not present in this sample\u2019s predictions.', true);
        return;
      }
      status(''); state.gene = v; hideGeneList(); render(true);
    }
    els.geneSearch.addEventListener('change', commitGene);
    els.geneSearch.addEventListener('keydown', function (e) {
      if (e.key === 'Enter') { e.preventDefault(); commitGene(); }
    });
    els.geneSearch.addEventListener('input', function () {
      if (geneListVisible()) renderGeneList(els.geneSearch.value);
      clearTimeout(timer);
      timer = setTimeout(function () {
        var v = els.geneSearch.value.trim();
        if (state.meta && inSample(v) && v !== state.gene) commitGene();
      }, 220);
    });

    document.addEventListener('keydown', function (e) {
      if (document.activeElement === els.geneSearch) return;
      if (!state.meta) return;
      if (e.key === 'ArrowRight' || e.key === 'n') step(1);
      else if (e.key === 'ArrowLeft' || e.key === 'p') step(-1);
      else if (e.key === 'f') { setWorld(); drawAll(); }
    });

    window.addEventListener('resize', function () { clampCamera(); drawAll(); });
    window.addEventListener('hashchange', function () {
      var h = readHash();
      if (h.sample && h.sample !== state.sample) selectSample(h.sample, h);
    });
  }

  function step(d) {
    var g = availGenes();
    var i = g.indexOf(state.gene);
    if (i < 0) i = 0;
    state.gene = g[(i + d + g.length) % g.length];
    els.geneSearch.value = state.gene;
    render(true);
  }

  function init() {
    ids.forEach(function (id) { els[id] = $(id); });
    attachPanZoom(els.heCanvas);
    attachPanZoom(els.predCanvas);
    els.cmapSel.value = state.cmap;
    wire();
    // the LUAD gene list used for the picker; the sample's own genes remain the
    // source of truth for what can actually be rendered
    fetch('/static/LUAD_gene.csv').then(function (r) {
      if (!r.ok) throw new Error('HTTP ' + r.status);
      return r.text();
    }).then(function (txt) {
      state.catalog = parseGeneCSV(txt);
    }).catch(function () {
      state.catalog = [];
    });
    api('/api/health').then(function (h) {
      els.backendInfo.textContent = (h.backend === 'hf' ? 'Hugging Face (on demand)'
        : 'local files')
        + (h.he_images && h.he_images.indexed ? ' \u00b7 ' + h.he_images.indexed + ' H&E' : '');
    }).catch(function () {});
    loadSamples();
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();





