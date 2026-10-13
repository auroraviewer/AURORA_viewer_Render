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
    celltype: null,
    kind: 'gene',          // 'gene' | 'celltype'
    cell: 14,
    cmap: 'turbo',
    opacity: 1.0,
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
  var ids = ['sampleSel', 'kindSel', 'geneSearch', 'geneList', 'celltypeSel', 'celltypeField',
             'geneField', 'resSeg', 'cmapSel', 'opacity', 'renderBtn', 'stage',
             'heCanvas', 'predCanvas', 'spinner', 'status',
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
  function busy(on) { els.spinner.style.display = on ? 'block' : 'none'; }
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
    ctx.globalAlpha = isHE ? 1 : state.opacity;

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
        ? h.sample : d.samples[0].sample_id;
      els.sampleSel.value = want;
      return selectSample(want, h);
    }).catch(function (e) { status('Cannot list samples: ' + e.message, true); });
  }

  function selectSample(id, restore) {
    state.sample = id;
    busy(true); status('Reading ' + id + '\u2026');
    return api('/api/samples/' + encodeURIComponent(id) + '/meta').then(function (m) {
      state.meta = m;

      // gene + cell-type pickers
      els.geneList.innerHTML = '';
      var frag = document.createDocumentFragment();
      m.genes.forEach(function (g) {
        var o = document.createElement('option'); o.value = g; frag.appendChild(o);
      });
      els.geneList.appendChild(frag);

      els.celltypeSel.innerHTML = '';
      m.cell_types.forEach(function (c) {
        var o = document.createElement('option'); o.value = c; o.textContent = c;
        els.celltypeSel.appendChild(o);
      });

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

      var gene = (restore && restore.gene && m.genes.indexOf(restore.gene) >= 0)
        ? restore.gene : m.genes[0];
      els.geneSearch.value = gene;
      state.gene = gene;
      state.celltype = m.cell_types[0];
      els.celltypeSel.value = state.celltype;
      if (restore && restore.kind) setKind(restore.kind, true);
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
    setWorld();
    render(true);
  }

  function syncResButtons() {
    Array.prototype.forEach.call(els.resSeg.children, function (b, i) {
      var L = state.meta.layers[i];
      b.setAttribute('aria-pressed', String(L && L.cell_px === state.cell));
    });
  }

  function setKind(kind, silent) {
    state.kind = kind;
    els.kindSel.value = kind;
    els.geneField.style.display = kind === 'gene' ? '' : 'none';
    els.celltypeField.style.display = kind === 'celltype' ? '' : 'none';
    if (!silent) render(false);
  }

  function render(keepView) {
    if (!state.sample) return;
    var name = state.kind === 'gene' ? state.gene : state.celltype;
    if (!name) return;
    var url = (state.kind === 'gene'
      ? '/api/gene.png?sample=' + encodeURIComponent(state.sample)
        + '&gene=' + encodeURIComponent(name)
      : '/api/celltype.png?sample=' + encodeURIComponent(state.sample)
        + '&celltype=' + encodeURIComponent(name))
      + '&cell=' + state.cell + '&width=2080';

    busy(true);
    status('Rendering ' + name + ' at ' + (currentLayer() ? currentLayer().microns : '?') + ' \u00b5m\u2026');
    var t0 = performance.now();

    // ask where this view sits in image-pixel space, in parallel with the image
    var statsUrl = '/api/view/stats?sample=' + encodeURIComponent(state.sample)
      + '&name=' + encodeURIComponent(name)
      + '&kind=' + state.kind + '&cell=' + state.cell + '&width=2080';

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
    var name = state.kind === 'gene' ? state.gene : state.celltype;
    var ri = state.renderInfo;
    if (state.kind === 'gene') {
      els.legendTitle.textContent = 'Predicted expression';
    } else {
      els.legendTitle.textContent = 'Cell-type proportion';
    }
    var L = currentLayer();
    var extra = '';
    if (state.kind === 'celltype' && ri && ri.raw_lo != null && ri.raw_hi != null) {
      // colour is the raw probability, stretched over the range this cell type
      // actually occupies; quote those real numbers so the colours are readable
      extra = ' \u00b7 ' + (100 * ri.raw_lo).toFixed(1) + '% to '
            + (100 * ri.raw_hi).toFixed(1) + '% of cells';
    } else if (state.kind === 'celltype') {
      extra = ' \u00b7 fraction of cells';
    }
    els.legendCaption.textContent = name
      + (L ? ' \u00b7 ' + L.microns + ' \u00b5m' : '')
      + extra;
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
      + '&kind=' + state.kind
      + (state.kind === 'gene' ? '&gene=' + encodeURIComponent(state.gene)
                               : '&celltype=' + encodeURIComponent(state.celltype))
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
    if (out.kind) state.kind = out.kind;
    return out;
  }

  // ---------------------------------------------------- wiring
  function wire() {
    els.sampleSel.addEventListener('change', function () { selectSample(els.sampleSel.value, null); });
    els.kindSel.addEventListener('change', function () { setKind(els.kindSel.value); });
    els.celltypeSel.addEventListener('change', function () {
      state.celltype = els.celltypeSel.value; render(true);
    });
    els.cmapSel.addEventListener('change', function () {
      state.cmap = els.cmapSel.value;
      if (state.bitmap) { buildValues(state.bitmap); drawAll(); }
      updateLegend(); updateUrl();
    });
    els.opacity.addEventListener('input', function () {
      state.opacity = parseInt(els.opacity.value, 10) / 100;
      drawAll();
    });
    els.renderBtn.addEventListener('click', function () { setWorld(); render(false); });

    var timer = null;
    function commitGene() {
      var v = els.geneSearch.value.trim();
      if (!v || !state.meta) return;
      if (state.meta.genes.indexOf(v) < 0) {
        status('Unknown gene "' + v + '" - pick one from the list.', true);
        return;
      }
      status(''); state.gene = v; render(true);
    }
    els.geneSearch.addEventListener('change', commitGene);
    els.geneSearch.addEventListener('keydown', function (e) {
      if (e.key === 'Enter') { e.preventDefault(); commitGene(); }
    });
    els.geneSearch.addEventListener('input', function () {
      clearTimeout(timer);
      timer = setTimeout(function () {
        var v = els.geneSearch.value.trim();
        if (state.meta && state.meta.genes.indexOf(v) >= 0 && v !== state.gene) commitGene();
      }, 220);
    });

    document.addEventListener('keydown', function (e) {
      if (document.activeElement === els.geneSearch) return;
      if (!state.meta || state.kind !== 'gene') return;
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
    var g = state.meta.genes;
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
    els.opacity.value = String(Math.round(state.opacity * 100));
    wire();
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





