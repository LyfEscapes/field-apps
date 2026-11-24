/*!
 * SiteLock 1.0 · LyfEscapes field camera kit
 * Return to the exact same spot, at the exact same angle, and measure in real millimetres.
 * Each app ships its own copy, so every app runs on its own. Apps on the same domain share
 * one on-device database ("lyf-field") and the Field Record format, so they also work together.
 * Marker detection: js-aruco2 (MIT), loaded from vendor/.
 */
(function (G) {
'use strict';
const SL = { version: '1.2.0' };
const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
const D2R = Math.PI / 180, R2D = 180 / Math.PI;
SL.clamp = clamp;

/* ---------------------------------------------------------------- homography */
function solve(A, b) {
  const n = b.length, M = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) {
    let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    [M[c], M[p]] = [M[p], M[c]];
    const v = M[c][c]; if (Math.abs(v) < 1e-12) return null;
    for (let r = 0; r < n; r++) { if (r === c) continue; const f = M[r][c] / v; if (!f) continue; for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k]; }
  }
  return M.map((r, i) => r[n] / r[i]);
}
SL.H = {
  mul(A, B) { const C = new Array(9); for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) C[i * 3 + j] = A[i * 3] * B[j] + A[i * 3 + 1] * B[3 + j] + A[i * 3 + 2] * B[6 + j]; return C; },
  inv(m) {
    const [a, b, c, d, e, f, g, h, i] = m, A = e * i - f * h, B = -(d * i - f * g), C = d * h - e * g, det = a * A + b * B + c * C;
    if (Math.abs(det) < 1e-15) return null;
    return [A / det, -(b * i - c * h) / det, (b * f - c * e) / det, B / det, (a * i - c * g) / det, -(a * f - c * d) / det, C / det, -(a * h - b * g) / det, (a * e - b * d) / det];
  },
  apply(H, p) { const w = H[6] * p[0] + H[7] * p[1] + H[8]; return [(H[0] * p[0] + H[1] * p[1] + H[2]) / w, (H[3] * p[0] + H[4] * p[1] + H[5]) / w]; },
  fit(src, dst) { // least squares, 4+ pairs of [x,y]
    if (src.length < 4) return null;
    const norm = P => { let mx = 0, my = 0; P.forEach(p => { mx += p[0]; my += p[1]; }); mx /= P.length; my /= P.length; let d = 0; P.forEach(p => d += Math.hypot(p[0] - mx, p[1] - my)); d /= P.length; const s = Math.SQRT2 / (d || 1); return { T: [s, 0, -s * mx, 0, s, -s * my, 0, 0, 1], P: P.map(p => [(p[0] - mx) * s, (p[1] - my) * s]) }; };
    const a = norm(src), b = norm(dst), M = Array.from({ length: 8 }, () => new Array(8).fill(0)), v = new Array(8).fill(0);
    for (let i = 0; i < src.length; i++) {
      const [x, y] = a.P[i], [u, w] = b.P[i], r1 = [x, y, 1, 0, 0, 0, -u * x, -u * y], r2 = [0, 0, 0, x, y, 1, -w * x, -w * y];
      for (let j = 0; j < 8; j++) { for (let k = 0; k < 8; k++) M[j][k] += r1[j] * r1[k] + r2[j] * r2[k]; v[j] += r1[j] * u + r2[j] * w; }
    }
    const h = solve(M, v); if (!h) return null;
    const Bi = SL.H.inv(b.T); return Bi && SL.H.mul(Bi, SL.H.mul([...h, 1], a.T));
  }
};

/* ------------------------------------------- similarity transform (ref -> current) */
// cur = C + s·R(rot)·(ref − C) + d, with d as a fraction of image width on both axes
SL.T = {
  id: () => ({ dx: 0, dy: 0, rot: 0, scale: 1 }),
  apply(t, p, w, h) { const cx = w / 2, cy = h / 2, c = Math.cos(t.rot * D2R) * t.scale, s = Math.sin(t.rot * D2R) * t.scale, x = p[0] - cx, y = p[1] - cy; return [cx + c * x - s * y + t.dx * w, cy + s * x + c * y + t.dy * w]; },
  invert(t, p, w, h) { const cx = w / 2, cy = h / 2, x = p[0] - cx - t.dx * w, y = p[1] - cy - t.dy * w, c = Math.cos(-t.rot * D2R) / t.scale, s = Math.sin(-t.rot * D2R) / t.scale; return [cx + c * x - s * y, cy + s * x + c * y]; },
  toH(t, w, h) { const cx = w / 2, cy = h / 2, c = Math.cos(t.rot * D2R) * t.scale, s = Math.sin(t.rot * D2R) * t.scale; return [c, -s, cx - c * cx + s * cy + t.dx * w, s, c, cy - s * cx - c * cy + t.dy * w, 0, 0, 1]; },
  fromH(H, w, h) { const C = [w / 2, h / 2], a = SL.H.apply(H, C), b = SL.H.apply(H, [C[0] + w / 4, C[1]]), vx = b[0] - a[0], vy = b[1] - a[1]; return { dx: (a[0] - C[0]) / w, dy: (a[1] - C[1]) / w, rot: Math.atan2(vy, vx) * R2D, scale: Math.hypot(vx, vy) / (w / 4) }; },
  // draw an image taken with transform t so it lines up with the reference frame (w × h)
  drawAligned(ctx, img, t, w, h) { ctx.save(); ctx.translate(w / 2, h / 2); ctx.rotate(-t.rot * D2R); ctx.scale(1 / t.scale, 1 / t.scale); ctx.translate(-w / 2 - t.dx * w, -h / 2 - t.dy * w); ctx.drawImage(img, 0, 0, w, h); ctx.restore(); }
};

// warp the current picture into the reference frame using H (ref px -> cur px), bilinear
SL.warp = (img, H, W, Hh, maxW = 1200) => {
  const k = Math.min(1, maxW / W), w = Math.round(W * k), h = Math.round(Hh * k), [sw, sh] = dims(img);
  const sc = document.createElement('canvas'); sc.width = sw; sc.height = sh; const sg = sc.getContext('2d', { willReadFrequently: true }); sg.drawImage(img, 0, 0); const S = sg.getImageData(0, 0, sw, sh).data;
  const out = document.createElement('canvas'); out.width = w; out.height = h; const og = out.getContext('2d'), O = og.createImageData(w, h), od = O.data;
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) {
    const X = x / k, Y = y / k, ww = H[6] * X + H[7] * Y + H[8], sx = (H[0] * X + H[1] * Y + H[2]) / ww, sy = (H[3] * X + H[4] * Y + H[5]) / ww, o = (y * w + x) * 4;
    if (sx < 0 || sy < 0 || sx >= sw - 1 || sy >= sh - 1) { od[o + 3] = 255; continue; }
    const x0 = sx | 0, y0 = sy | 0, fx = sx - x0, fy = sy - y0, i = (y0 * sw + x0) * 4, j = i + sw * 4;
    for (let c = 0; c < 3; c++) od[o + c] = S[i + c] * (1 - fx) * (1 - fy) + S[i + 4 + c] * fx * (1 - fy) + S[j + c] * (1 - fx) * fy + S[j + 4 + c] * fx * fy;
    od[o + 3] = 255;
  }
  og.putImageData(O, 0, 0); return out;
};

/* ------------------------------------------------------------------- FFT core */
function fft1(re, im, n, inv) {
  for (let i = 1, j = 0; i < n; i++) { let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit; if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; } }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = 2 * Math.PI / len * (inv ? 1 : -1), wr = Math.cos(ang), wi = Math.sin(ang), half = len >> 1;
    for (let i = 0; i < n; i += len) {
      let cr = 1, ci = 0;
      for (let k = 0; k < half; k++) {
        const a = i + k, b = a + half, xr = re[b] * cr - im[b] * ci, xi = re[b] * ci + im[b] * cr;
        re[b] = re[a] - xr; im[b] = im[a] - xi; re[a] += xr; im[a] += xi;
        const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
      }
    }
  }
  if (inv) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
}
function fft2(re, im, N, inv) {
  const rr = new Float64Array(N), ri = new Float64Array(N);
  for (let y = 0; y < N; y++) { const o = y * N; for (let x = 0; x < N; x++) { rr[x] = re[o + x]; ri[x] = im[o + x]; } fft1(rr, ri, N, inv); for (let x = 0; x < N; x++) { re[o + x] = rr[x]; im[o + x] = ri[x]; } }
  for (let x = 0; x < N; x++) { for (let y = 0; y < N; y++) { rr[y] = re[y * N + x]; ri[y] = im[y * N + x]; } fft1(rr, ri, N, inv); for (let y = 0; y < N; y++) { re[y * N + x] = rr[y]; im[y * N + x] = ri[y]; } }
}

/* ------------------------------------------------------------- alignment engine */
// Gaussian low-pass on the cross-power spectrum: smooths the match surface so small turn/scale errors don't swamp it
const LP = {};
function lowpass(N) { if (LP[N]) return LP[N]; const w = new Float64Array(N * N), sg = N * .16; let sum = 0; for (let y = 0; y < N; y++) { const fy = y < N / 2 ? y : y - N; for (let x = 0; x < N; x++) { const fx = x < N / 2 ? x : x - N, v = Math.exp(-(fx * fx + fy * fy) / (2 * sg * sg)); w[y * N + x] = v; sum += v; } } w.mean = sum / (N * N); return LP[N] = w; }
const scratch = {};
function canvasOf(key, w, h) { let c = scratch[key]; if (!c) { c = scratch[key] = document.createElement('canvas'); } if (c.width !== w) c.width = w; if (c.height !== h) c.height = h; return c; }
const dims = s => [s.videoWidth || s.naturalWidth || s.width, s.videoHeight || s.naturalHeight || s.height];
SL.Align = {
  // turn a picture into a windowed edge map ready for phase correlation
  prep(src, N = 128) {
    const [sw, sh] = dims(src), k = N / Math.max(sw, sh), w = Math.round(sw * k), h = Math.round(sh * k), ox = (N - w) >> 1, oy = (N - h) >> 1;
    const c = canvasOf('prep' + N, N, N), g = c.getContext('2d', { willReadFrequently: true });
    g.fillStyle = '#808080'; g.fillRect(0, 0, N, N); g.drawImage(src, ox, oy, w, h);
    const px = g.getImageData(0, 0, N, N).data, L = new Float64Array(N * N), E = new Float64Array(N * N);
    for (let i = 0; i < N * N; i++) L[i] = px[i * 4] * .299 + px[i * 4 + 1] * .587 + px[i * 4 + 2] * .114;
    let sum = 0, n = 0;
    for (let y = oy + 1; y < oy + h - 1; y++) for (let x = ox + 1; x < ox + w - 1; x++) { const i = y * N + x, v = Math.abs(L[i + 1] - L[i - 1]) + Math.abs(L[i + N] - L[i - N]); E[i] = v; sum += v; n++; }
    const mean = n ? sum / n : 0;
    for (let y = oy + 1; y < oy + h - 1; y++) { const wy = .5 - .5 * Math.cos(2 * Math.PI * (y - oy) / (h - 1)); for (let x = ox + 1; x < ox + w - 1; x++) { const wx = .5 - .5 * Math.cos(2 * Math.PI * (x - ox) / (w - 1)), i = y * N + x; E[i] = (E[i] - mean) * wx * wy; } }
    return { N, E, k, sw, sh, F: null };
  },
  spectrum(P) { if (!P.F) { const re = Float64Array.from(P.E), im = new Float64Array(P.E.length); fft2(re, im, P.N, false); P.F = [re, im]; } return P.F; },
  warp(P, rot, scale) { // sample P at C + s·R·(p − C)
    const N = P.N, out = new Float64Array(N * N), c = Math.cos(rot * D2R) * scale, s = Math.sin(rot * D2R) * scale, h = N / 2, E = P.E;
    for (let y = 0; y < N; y++) { const yy = y - h; for (let x = 0; x < N; x++) { const xx = x - h, sx = h + c * xx - s * yy, sy = h + s * xx + c * yy, x0 = sx | 0, y0 = sy | 0; if (x0 < 0 || y0 < 0 || x0 >= N - 1 || y0 >= N - 1) continue; const fx = sx - x0, fy = sy - y0, i = y0 * N + x0; out[y * N + x] = E[i] * (1 - fx) * (1 - fy) + E[i + 1] * fx * (1 - fy) + E[i + N] * (1 - fx) * fy + E[i + N + 1] * fx * fy; } }
    return out;
  },
  correlate(refP, curP, rot, scale) {
    const N = refP.N, [Fr, Fi] = SL.Align.spectrum(refP), re = SL.Align.warp(curP, rot, scale), im = new Float64Array(N * N);
    fft2(re, im, N, false);
    const Wt = lowpass(N);
    for (let i = 0; i < N * N; i++) { const a = re[i], b = im[i], c = Fr[i], d = -Fi[i], r = a * c - b * d, j = a * d + b * c, m = (Math.hypot(r, j) + 1e-9) / Wt[i]; re[i] = r / m; im[i] = j / m; }
    fft2(re, im, N, true);
    for (let i = 0; i < N * N; i++) re[i] /= Wt.mean;
    let best = -1, bi = 0; for (let i = 0; i < N * N; i++) if (re[i] > best) { best = re[i]; bi = i; }
    const px = bi % N, py = (bi / N) | 0, at = (x, y) => re[((y + N) % N) * N + ((x + N) % N)];
    const sub = (m, c, p) => { const d = m - 2 * c + p; return d ? clamp(.5 * (m - p) / d, -.5, .5) : 0; };
    let ex = px + sub(at(px - 1, py), best, at(px + 1, py)), ey = py + sub(at(px, py - 1), best, at(px, py + 1));
    if (ex > N / 2) ex -= N; if (ey > N / 2) ey -= N;
    return { peak: best, ex, ey };
  },
  // find how the current picture is shifted, turned and scaled against the reference
  register(refP, curP, opt = {}) {
    const near = opt.near, rots = [], scs = [];
    if (near) { for (let r = -2; r <= 2; r += 1) rots.push(near.rot + r); for (let s = -2; s <= 2; s++) scs.push(near.scale * (1 + s * .02)); }
    else { for (let r = -15; r <= 15; r += 3) rots.push(r); for (const s of [.8, .87, .94, 1, 1.07, 1.15, 1.25]) scs.push(s); }
    let best = null;
    for (const r of rots) for (const s of scs) { const m = SL.Align.correlate(refP, curP, r, s); if (!best || m.peak > best.peak) best = { ...m, rot: r, scale: s }; }
    if (opt.refine !== false) { // local refinement
      const b0 = best;
      for (const dr of [-1, -.5, .5, 1]) { const m = SL.Align.correlate(refP, curP, b0.rot + dr, b0.scale); if (m.peak > best.peak) best = { ...m, rot: b0.rot + dr, scale: b0.scale }; }
      const b1 = best;
      for (const ds of [-.015, -.0075, .0075, .015]) { const m = SL.Align.correlate(refP, curP, b1.rot, b1.scale * (1 + ds)); if (m.peak > best.peak) best = { ...m, rot: b1.rot, scale: b1.scale * (1 + ds) }; }
    }
    // warped(p) = ref(p − e) → d = s·R·e ; convert from prep pixels to fraction of image width
    const c = Math.cos(best.rot * D2R) * best.scale, s = Math.sin(best.rot * D2R) * best.scale, dxp = c * best.ex - s * best.ey, dyp = s * best.ex + c * best.ey;
    const wpx = curP.sw * curP.k;
    return { dx: dxp / wpx, dy: dyp / wpx, rot: best.rot, scale: best.scale, peak: best.peak, q: clamp((best.peak - .2) / .35, 0, 1) };
  },
  // full-quality match for a saved photo: coarse at 128, then refined at 256
  registerStill(refImg, curImg) {
    const a = SL.Align.prep(refImg, 128), b = SL.Align.prep(curImg, 128), coarse = SL.Align.register(a, b);
    const A = SL.Align.prep(refImg, 256), B = SL.Align.prep(curImg, 256);
    const fine = SL.Align.register(A, B, { near: { rot: coarse.rot, scale: coarse.scale } });
    return fine.peak >= coarse.peak * .8 ? fine : coarse;
  },
  // exact match from the reference card seen in both pictures: ref px -> card mm -> cur px
  fromCards(refCardH, curCardH) { const Hi = SL.H.inv(curCardH); return Hi && SL.H.mul(Hi, refCardH); },
  // best available alignment of a saved photo: card if both have it (checked against the wall texture), else the image match
  still(refImg, curImg, refCard, curCard) {
    const [w, h] = dims(refImg), img = SL.Align.registerStill(refImg, curImg), Himg = SL.T.toH(img, w, h);
    if (refCard && curCard) {
      const Hc = SL.Align.fromCards(refCard.H, curCard.H);
      if (Hc) { const tc = SL.T.fromH(Hc, w, h), moved = img.q > .5 && (Math.hypot(tc.dx - img.dx, tc.dy - img.dy) > .03 || Math.abs(tc.scale - img.scale) > .05);
        if (!moved) return { H: Hc, t: { ...tc, q: Math.max(img.q, .9) }, method: 'card', q: Math.max(img.q, .9) };
        return { H: Himg, t: img, method: 'image', q: img.q, cardMoved: true }; }
    }
    return { H: Himg, t: img, method: 'image', q: img.q };
  },
  // plain-language advice; sim=true phrases it as drag directions
  advice(t, sim) {
    const out = [], sh = Math.hypot(t.dx, t.dy);
    if (Math.abs(t.scale - 1) > .035) out.push({ k: 'z', v: Math.abs(t.scale - 1) * 2, text: sim ? (t.scale > 1 ? 'Pinch to shrink' : 'Spread to enlarge') : (t.scale > 1 ? 'Step back a little' : 'Move a little closer') });
    if (sh > .02) { const h = Math.abs(t.dx) > Math.abs(t.dy) * 1.6 ? (t.dx > 0 ? 'right' : 'left') : Math.abs(t.dy) > Math.abs(t.dx) * 1.6 ? (t.dy > 0 ? 'down' : 'up') : `${t.dy > 0 ? 'down' : 'up'} and ${t.dx > 0 ? 'right' : 'left'}`; const opp = { right: 'left', left: 'right', up: 'down', down: 'up' }; out.push({ k: 'm', v: sh * 3, text: sim ? 'Drag ' + h.replace(/right|left|up|down/g, m => opp[m]) : 'Move ' + h }); }
    if (Math.abs(t.rot) > 1.5) out.push({ k: 'r', v: Math.abs(t.rot) / 10, text: sim ? (t.rot > 0 ? 'Rotate left' : 'Rotate right') : (t.rot > 0 ? 'Turn the phone clockwise' : 'Turn the phone counter-clockwise') });
    out.sort((a, b) => b.v - a.v);
    return out;
  },
  locked: t => Math.hypot(t.dx, t.dy) < .02 && Math.abs(t.rot) < 1.5 && Math.abs(t.scale - 1) < .035 && t.q > .35
};

/* ------------------------------------------------------ reference card & spot tags */
// Card: credit-card size, four ARUCO_MIP_36h12 markers. Card serial k uses marker ids 4k..4k+3 (serials 0–39).
// Field board: US-letter landscape with 50 mm markers for culverts, headwalls and ditches (serials 40–49).
const CARD = { kind: 'card', W: 85.6, H: 54, M: 12, at: [[4, 4], [69.6, 4], [69.6, 38], [4, 38]] };
const BOARD = { kind: 'board', W: 279.4, H: 215.9, M: 50, at: [[10, 10], [219.4, 10], [219.4, 155.9], [10, 155.9]] };
const specOf = serial => serial >= 40 ? BOARD : CARD;
let detector = null;
const getDetector = () => detector || (G.AR ? (detector = new G.AR.Detector({ dictionaryName: 'ARUCO_MIP_36h12', maxHammingDistance: 4 })) : null);
const markerSVG = (id, x, y, size) => { const d = G.AR && new G.AR.Dictionary('ARUCO_MIP_36h12'); if (!d) return `<rect x="${x}" y="${y}" width="${size}" height="${size}" fill="#000"/>`; const q = size * 10 / 8, o = size / 8; return d.generateSVG(id).replace('<svg ', `<svg x="${(x - o).toFixed(3)}" y="${(y - o).toFixed(3)}" width="${q.toFixed(3)}" height="${q.toFixed(3)}" shape-rendering="crispEdges" `); };
SL.Card = {
  spec: CARD, board: BOARD, specOf,
  boardSVG(serial = 40, brand = 'SITELOCK') {
    let t = ''; for (let i = 0; i <= 200; i++) { const h = i % 10 === 0 ? 7 : i % 5 === 0 ? 5 : 3; t += `<rect x="${(39.7 + i - .2).toFixed(2)}" y="${(30 - h).toFixed(2)}" width=".4" height="${h}"/>`; }
    let u = ''; for (let i = 0; i <= 64; i++) { const x = 39.7 + i * 25.4 / 8, h = i % 8 === 0 ? 7 : i % 4 === 0 ? 5 : 3; u += `<rect x="${(x - .2).toFixed(2)}" y="185" width=".4" height="${h}"/>`; if (i % 8 === 0) u += `<text x="${x.toFixed(2)}" y="199" font-size="5" text-anchor="middle">${i / 8}</text>`; }
    const marks = BOARD.at.map((p, j) => markerSVG(serial * 4 + j, p[0], p[1], BOARD.M)).join('');
    return `<svg xmlns="http://www.w3.org/2000/svg" width="279.4mm" height="215.9mm" viewBox="0 0 279.4 215.9" font-family="Arial,Helvetica,sans-serif"><rect x=".5" y=".5" width="278.4" height="214.9" fill="#fff" stroke="#c9d2df"/>${marks}<g fill="#0b1424"><rect x="39.7" y="30" width="200" height=".8"/>${t}<text x="39.7" y="38" font-size="5">0</text><text x="239.7" y="38" font-size="5" text-anchor="end">200 mm</text><text x="139.7" y="95" font-size="16" font-weight="700" text-anchor="middle">${brand} FIELD BOARD</text><text x="139.7" y="112" font-size="8" text-anchor="middle">Board ${String(serial).padStart(3, '0')} · prop it flat against the face you are measuring</text><text x="139.7" y="124" font-size="6" text-anchor="middle">Print on US letter at 100%, check the 200 mm bar, laminate or mount on plywood</text><rect x="39.7" y="184" width="203.2" height=".8"/>${u}<text x="242.9" y="181" font-size="5" text-anchor="end">inches</text></g></svg>`;
  },
  ready: () => !!getDetector(),
  svg(serial = 0, brand = 'SITELOCK') {
    if (serial >= 40) return SL.Card.boardSVG(serial, brand);
    let t = ''; for (let i = 0; i <= 50; i++) { const h = i % 10 === 0 ? 3.2 : i % 5 === 0 ? 2.2 : 1.3; t += `<rect x="${(18 + i - .1).toFixed(2)}" y="${(10 - h).toFixed(2)}" width=".2" height="${h}"/>`; }
    const cmp = [.1, .2, .3, .5, 1, 2, 3].map((w, i) => `<rect x="${20 + i * 7}" y="23" width="${w}" height="8"/><text x="${20 + i * 7}" y="34.5" font-size="2.1">${w}</text>`).join('');
    const marks = CARD.at.map((p, j) => markerSVG(serial * 4 + j, p[0], p[1], CARD.M)).join('');
    return `<svg xmlns="http://www.w3.org/2000/svg" width="85.6mm" height="54mm" viewBox="0 0 85.6 54" font-family="Arial,Helvetica,sans-serif"><rect x=".2" y=".2" width="85.2" height="53.6" rx="3.2" fill="#fff" stroke="#c9d2df" stroke-width=".3"/>${marks}<g fill="#0b1424"><rect x="18" y="10" width="50" height=".45"/>${t}<text x="18" y="14.2" font-size="2.2">0</text><text x="68" y="14.2" font-size="2.2" text-anchor="end">50 mm</text>${cmp}<text x="68" y="34.5" font-size="2.1" text-anchor="end">mm</text><text x="42.8" y="44" font-size="3.6" font-weight="700" text-anchor="middle">${brand}</text><text x="42.8" y="48.5" font-size="2.3" text-anchor="middle">Card ${String(serial).padStart(3, '0')} · keep flat beside the target</text></g></svg>`;
  },
  tagSVG(id, label) { // single-marker spot sticker, ids 200–249, 20 mm marker
    return `<svg xmlns="http://www.w3.org/2000/svg" width="30mm" height="38mm" viewBox="0 0 30 38" font-family="Arial,Helvetica,sans-serif"><rect x=".2" y=".2" width="29.6" height="37.6" rx="2.5" fill="#fff" stroke="#c9d2df" stroke-width=".3"/>${markerSVG(id, 5, 5, 20)}<text x="15" y="31" font-size="3.4" font-weight="700" text-anchor="middle" fill="#0b1424">SPOT ${id}</text><text x="15" y="35" font-size="2.2" text-anchor="middle" fill="#0b1424">${label || 'Leave in place'}</text></svg>`;
  },
  // find markers; returns card pose (image px -> card mm) and/or spot tags
  detect(src, maxSide = 1000) {
    const det = getDetector(); if (!det) return null;
    const [sw, sh] = dims(src), k = Math.min(1, maxSide / Math.max(sw, sh)), w = Math.round(sw * k), h = Math.round(sh * k);
    const c = canvasOf('det', w, h), g = c.getContext('2d', { willReadFrequently: true }); g.drawImage(src, 0, 0, w, h);
    let ms; try { ms = det.detect(g.getImageData(0, 0, w, h)); } catch (e) { return null; }
    const up = p => [p.x / k, p.y / k], tags = [], groups = {};
    ms.forEach(m => { if (m.id >= 200 && m.id < 250) tags.push({ id: m.id, corners: m.corners.map(up) }); else if (m.id < 248) (groups[m.id >> 2] = groups[m.id >> 2] || []).push(m); });
    const out = { w: sw, h: sh, tags, card: null };
    const serial = Object.keys(groups).sort((a, b) => groups[b].length - groups[a].length)[0];
    if (serial != null) {
      const img = [], mm = [];
      const SP = specOf(+serial);
      groups[serial].forEach(m => { const [x, y] = SP.at[m.id & 3], M = SP.M, ref = [[x, y], [x + M, y], [x + M, y + M], [x, y + M]]; m.corners.forEach((p, i) => { img.push(up(p)); mm.push(ref[i]); }); });
      const H = SL.H.fit(img, mm), Hi = H && SL.H.inv(H);
      if (H && Hi) {
        const quad = [[0, 0], [SP.W, 0], [SP.W, SP.H], [0, SP.H]].map(p => SL.H.apply(Hi, p));
        const ctr = SL.H.apply(Hi, [SP.W / 2, SP.H / 2]), a = SL.H.apply(H, ctr), b = SL.H.apply(H, [ctr[0] + 1, ctr[1]]);
        out.card = { kind: SP.kind, serial: +serial, markers: groups[serial].length, H, Hi, quad, mmPerPx: Math.hypot(b[0] - a[0], b[1] - a[1]), tilt: SL.Card.tilt(Hi, sw, sh) };
      }
    }
    return out;
  },
  // angle between the camera's view and the card plane, from the homography (focal length estimated)
  tilt(Hi, w, h, fRatio = .78) {
    const f = fRatio * Math.max(w, h), cx = w / 2, cy = h / 2;
    const col = j => [(Hi[j] - cx * Hi[6 + j]) / f, (Hi[3 + j] - cy * Hi[6 + j]) / f, Hi[6 + j]];
    const r1 = col(0), r2 = col(1), n1 = Math.hypot(...r1), n2 = Math.hypot(...r2);
    const a = r1.map(v => v / n1), b = r2.map(v => v / n2), r3 = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    return Math.acos(clamp(Math.abs(r3[2]) / Math.hypot(...r3), 0, 1)) * R2D;
  },
  mm(H, p, q) { const a = SL.H.apply(H, p), b = SL.H.apply(H, q); return Math.hypot(a[0] - b[0], a[1] - b[1]); }
};
// lengths for people: mm, or inches (feet and inches past two feet)
SL.len = (mm, units = 'mm') => { if (mm == null || isNaN(mm)) return '–'; if (units !== 'in') return (mm < 1 ? mm.toFixed(2) : mm < 100 ? mm.toFixed(1) : Math.round(mm)) + ' mm'; const i = mm / 25.4; if (i >= 24) { const f = Math.floor(i / 12), r = Math.round(i - f * 12); return r === 12 ? `${f + 1} ft` : `${f} ft ${r} in`; } return (i < 1 ? i.toFixed(2) : i.toFixed(1)) + ' in'; };

/* ------------------------------------------------------------- GPS & sensors */
SL.units = 'imperial'; // 'metric' for m/km
SL.dist = m => SL.units === 'metric' ? (m < 1000 ? Math.round(m) + ' m' : (m / 1000).toFixed(1) + ' km') : (m < 1609 ? Math.round(m * 3.281) + ' ft' : (m / 1609.34).toFixed(1) + ' mi');
SL.Geo = {
  dist(a, b) { const R = 6371000, p1 = a.lat * D2R, p2 = b.lat * D2R, dp = (b.lat - a.lat) * D2R, dl = (b.lon - a.lon) * D2R, x = Math.sin(dp / 2) ** 2 + Math.cos(p1) * Math.cos(p2) * Math.sin(dl / 2) ** 2; return 2 * R * Math.atan2(Math.sqrt(x), Math.sqrt(1 - x)); },
  bearing(a, b) { const p1 = a.lat * D2R, p2 = b.lat * D2R, dl = (b.lon - a.lon) * D2R; return (Math.atan2(Math.sin(dl) * Math.cos(p2), Math.cos(p1) * Math.sin(p2) - Math.sin(p1) * Math.cos(p2) * Math.cos(dl)) * R2D + 360) % 360; },
  compass: d => ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'][Math.round(d / 45) % 8],
  fmt: g => g ? `${g.lat.toFixed(6)}, ${g.lon.toFixed(6)}${g.acc ? ` ±${SL.dist(g.acc)}` : ''}` : 'No GPS',
  watch(cb) {
    if (!navigator.geolocation) { cb({ err: 'GPS not available' }); return () => {}; }
    let id; try { id = navigator.geolocation.watchPosition(p => cb({ lat: p.coords.latitude, lon: p.coords.longitude, acc: p.coords.accuracy, alt: p.coords.altitude, t: p.timestamp }), e => cb({ err: e.code === 1 ? 'GPS blocked' : 'GPS unavailable' }), { enableHighAccuracy: true, maximumAge: 2000, timeout: 20000 }); } catch (e) { cb({ err: 'GPS blocked' }); return () => {}; }
    return () => { try { navigator.geolocation.clearWatch(id); } catch (e) {} };
  }
};
SL.Orient = {
  // call straight from a tap: iOS asks permission here
  start(cb) {
    let stop = () => {}, on = true;
    const attach = () => {
      if (!on) return;
      const evt = 'ondeviceorientationabsolute' in G ? 'deviceorientationabsolute' : 'deviceorientation';
      const h = e => { if (e.beta == null) return; const heading = e.webkitCompassHeading != null ? e.webkitCompassHeading : (e.absolute && e.alpha != null ? (360 - e.alpha) % 360 : null); cb({ pitch: e.beta, roll: e.gamma, heading }); };
      G.addEventListener(evt, h); stop = () => G.removeEventListener(evt, h);
    };
    try {
      if (G.DeviceOrientationEvent && typeof DeviceOrientationEvent.requestPermission === 'function') DeviceOrientationEvent.requestPermission().then(r => r === 'granted' ? attach() : cb({ err: 'Motion blocked' })).catch(() => cb({ err: 'Motion blocked' }));
      else if (G.DeviceOrientationEvent) attach(); else cb({ err: 'No motion sensor' });
    } catch (e) { cb({ err: 'Motion blocked' }); }
    return () => { on = false; stop(); };
  }
};

// phone as an inclinometer: lay it flat on a surface; returns the surface's angle from level
SL.Incline = {
  angle(o) { if (!o || o.pitch == null) return null; return Math.acos(clamp(Math.cos(o.pitch * D2R) * Math.cos((o.roll || 0) * D2R), -1, 1)) * R2D; },
  ratio(deg) { return deg >= 89.9 ? 0 : 1 / Math.tan(deg * D2R); } // horizontal : vertical
};

/* ---------------------------------------------------------------------- EXIF */
SL.Exif = {
  read(buf) {
    try {
      const v = new DataView(buf); if (v.getUint16(0) !== 0xFFD8) return {};
      let o = 2;
      while (o < v.byteLength - 4) {
        const m = v.getUint16(o), len = v.getUint16(o + 2);
        if (m === 0xFFE1 && v.getUint32(o + 4) === 0x45786966) return SL.Exif.tiff(v, o + 10);
        if ((m & 0xFF00) !== 0xFF00) break; o += 2 + len;
      }
    } catch (e) {}
    return {};
  },
  tiff(v, t) {
    const le = v.getUint16(t) === 0x4949, u16 = p => v.getUint16(p, le), u32 = p => v.getUint32(p, le);
    const ifd = p => { const n = u16(p), r = {}; for (let i = 0; i < n; i++) { const e = p + 2 + i * 12; r[u16(e)] = { type: u16(e + 2), count: u32(e + 4), off: e + 8 }; } return r; };
    const val = en => { if (!en) return null; const sz = { 1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 7: 1, 10: 8 }[en.type] || 1, at = sz * en.count > 4 ? t + u32(en.off) : en.off;
      if (en.type === 2) { let s = ''; for (let i = 0; i < en.count - 1; i++) s += String.fromCharCode(v.getUint8(at + i)); return s; }
      if (en.type === 3) return u16(at); if (en.type === 4) return u32(at);
      if (en.type === 5 || en.type === 10) { const a = []; for (let i = 0; i < en.count; i++) a.push(u32(at + i * 8) / (u32(at + i * 8 + 4) || 1)); return en.count === 1 ? a[0] : a; }
      return null; };
    const i0 = ifd(t + u32(t + 4)), out = { make: val(i0[0x10F]), model: val(i0[0x110]), orientation: val(i0[0x112]) };
    if (i0[0x8769]) { const ex = ifd(t + val(i0[0x8769])); out.date = val(ex[0x9003]); out.focal = val(ex[0x920A]); out.focal35 = val(ex[0xA405]); }
    if (i0[0x8825]) { const gp = ifd(t + val(i0[0x8825])), dms = a => a && a[0] + a[1] / 60 + a[2] / 3600, la = dms(val(gp[2])), lo = dms(val(gp[4])); if (la != null && lo != null) { out.lat = val(gp[1]) === 'S' ? -la : la; out.lon = val(gp[3]) === 'W' ? -lo : lo; out.alt = val(gp[6]); } }
    return out;
  }
};

/* --------------------------------------------------------------------- grid */
SL.Grid = {
  cols: 8, rows: 6,
  label(fx, fy, cols = 8, rows = 6) { const c = clamp(Math.floor(fx * cols), 0, cols - 1), r = clamp(Math.floor(fy * rows), 0, rows - 1); return String.fromCharCode(65 + c) + (r + 1); },
  draw(g, x, y, w, h, opt = {}) {
    const cols = opt.cols || 8, rows = opt.rows || 6, fs = Math.max(11, Math.min(15, w / 40));
    g.save(); g.strokeStyle = opt.color || 'rgba(255,255,255,.55)'; g.lineWidth = 1; g.beginPath();
    for (let i = 1; i < cols; i++) { const xx = Math.round(x + w * i / cols) + .5; g.moveTo(xx, y); g.lineTo(xx, y + h); }
    for (let i = 1; i < rows; i++) { const yy = Math.round(y + h * i / rows) + .5; g.moveTo(x, yy); g.lineTo(x + w, yy); }
    g.stroke(); g.font = `700 ${fs}px system-ui,sans-serif`; g.fillStyle = 'rgba(255,255,255,.85)'; g.shadowColor = 'rgba(0,0,0,.7)'; g.shadowBlur = 3;
    for (let i = 0; i < cols; i++) g.fillText(String.fromCharCode(65 + i), x + w * (i + .5) / cols - fs / 3, y + fs + 2);
    for (let i = 0; i < rows; i++) g.fillText(String(i + 1), x + 4, y + h * (i + .5) / rows + fs / 3);
    if (opt.pin) { const cx = x + opt.pin.x * w, cy = y + opt.pin.y * h, cw = w / cols, ch = h / rows, c = Math.floor(opt.pin.x * cols), r = Math.floor(opt.pin.y * rows);
      g.shadowBlur = 0; g.fillStyle = 'rgba(255,45,85,.18)'; g.fillRect(x + c * cw, y + r * ch, cw, ch); g.strokeStyle = '#ff2d55'; g.lineWidth = 2.5; g.strokeRect(x + c * cw, y + r * ch, cw, ch);
      g.beginPath(); g.arc(cx, cy, 9, 0, 7); g.moveTo(cx - 16, cy); g.lineTo(cx + 16, cy); g.moveTo(cx, cy - 16); g.lineTo(cx, cy + 16); g.stroke(); }
    g.restore();
  }
};

/* ------------------------------------------------- shared on-device store */
// One database per domain. Every LyfEscapes app writes here with its own app id,
// so apps installed side by side can see and link each other's records.
SL.Store = function (app) { this.app = app; this.mem = { spots: new Map(), caps: new Map(), kv: new Map() }; this.db = null; };
SL.Store.prototype.open = function () {
  return new Promise(res => {
    try {
      const rq = indexedDB.open('lyf-field', 1);
      rq.onupgradeneeded = () => { const db = rq.result; db.createObjectStore('spots', { keyPath: 'id' }).createIndex('app', 'app'); db.createObjectStore('caps', { keyPath: 'id' }).createIndex('spot', 'spot'); db.createObjectStore('kv', { keyPath: 'k' }); };
      rq.onsuccess = () => { this.db = rq.result; res(this); }; rq.onerror = () => res(this); rq.onblocked = () => res(this);
    } catch (e) { res(this); }
  });
};
SL.Store.prototype._tx = function (st, mode, fn) {
  return new Promise(res => {
    if (!this.db) return res(fn(null));
    try { const tx = this.db.transaction(st, mode), os = tx.objectStore(st), r = fn(os); if (r && 'onsuccess' in r) { r.onsuccess = () => res(r.result); r.onerror = () => res(null); } else tx.oncomplete = () => res(r); }
    catch (e) { res(fn(null)); }
  });
};
SL.Store.prototype.put = function (st, o) { this.mem[st].set(o.id || o.k, o); return this._tx(st, 'readwrite', os => os ? os.put(o) : o); };
SL.Store.prototype.get = function (st, id) { return this._tx(st, 'readonly', os => os ? os.get(id) : this.mem[st].get(id)); };
SL.Store.prototype.del = function (st, id) { this.mem[st].delete(id); return this._tx(st, 'readwrite', os => os ? os.delete(id) : null); };
SL.Store.prototype.all = function (st) { return this._tx(st, 'readonly', os => os ? os.getAll() : [...this.mem[st].values()]).then(a => a || []); };

/* ------------------------------------------------------------- Field Record */
SL.Record = {
  pack(app, spot, caps) { return JSON.stringify({ fieldRecord: 1, app, exported: new Date().toISOString(), spot, caps }, null, 1); },
  unpack(text) { const o = JSON.parse(text); if (!o || o.fieldRecord !== 1 || !o.spot) throw new Error('Not a Field Record'); return o; }
};

/* ------------------------------------------------------------ image helpers */
SL.img = src => new Promise((res, rej) => { const i = new Image(); i.onload = () => res(i); i.onerror = rej; i.src = src; });
SL.toJPEG = (src, max = 1600, q = .86) => { const [w, h] = dims(src), k = Math.min(1, max / Math.max(w, h)), c = document.createElement('canvas'); c.width = Math.round(w * k); c.height = Math.round(h * k); c.getContext('2d').drawImage(src, 0, 0, c.width, c.height); return c.toDataURL('image/jpeg', q); };
SL.edges = (src, max = 480) => { // bright outline of the reference for the ghost "edges" mode
  const [w, h] = dims(src), k = max / Math.max(w, h), W = Math.round(w * k), H = Math.round(h * k), c = document.createElement('canvas'); c.width = W; c.height = H;
  const g = c.getContext('2d', { willReadFrequently: true }); g.drawImage(src, 0, 0, W, H); const d = g.getImageData(0, 0, W, H), p = d.data, L = new Float32Array(W * H), E = new Float32Array(W * H);
  for (let i = 0; i < W * H; i++) L[i] = p[i * 4] * .3 + p[i * 4 + 1] * .59 + p[i * 4 + 2] * .11;
  for (let y = 1; y < H - 1; y++) for (let x = 1; x < W - 1; x++) { const i = y * W + x; E[i] = Math.abs(L[i + 1] - L[i - 1]) + Math.abs(L[i + W] - L[i - W]); }
  const s = Float32Array.from(E).sort(), th = s[Math.floor(s.length * .9)] || 1;
  for (let i = 0; i < W * H; i++) { const on = E[i] > th; p[i * 4] = 40; p[i * 4 + 1] = 240; p[i * 4 + 2] = 255; p[i * 4 + 3] = on ? 230 : 0; }
  g.putImageData(d, 0, 0); return c;
};

/* ----------------------------------------------------------------- compare */
SL.Compare = {
  mount(box, a, b, Hm, opt = {}) {
    const W = a.naturalWidth || a.width, H = a.naturalHeight || a.height, cv = document.createElement('canvas'); cv.width = W; cv.height = H; cv.className = 'sl-cmp';
    cv.style.cssText = 'width:100%;height:auto;display:block;border-radius:12px;touch-action:none;background:#111';
    box.innerHTML = ''; box.appendChild(cv);
    const g = cv.getContext('2d'), st = { mode: opt.mode || 'swipe', pos: .5, fade: .5, flick: false, grid: !!opt.grid, pin: opt.pin, timer: null };
    const bb = SL.warp(b, Hm || [1, 0, 0, 0, 1, 0, 0, 0, 1], W, H, W);
    let diff = null;
    const draw = () => {
      g.clearRect(0, 0, W, H);
      if (st.mode === 'swipe') { g.drawImage(a, 0, 0, W, H); g.save(); g.beginPath(); g.rect(W * st.pos, 0, W, H); g.clip(); g.drawImage(bb, 0, 0); g.restore(); g.fillStyle = '#fff'; g.fillRect(W * st.pos - 2, 0, 4, H); g.beginPath(); g.arc(W * st.pos, H / 2, Math.max(14, W / 50), 0, 7); g.fill(); }
      else if (st.mode === 'fade') { g.drawImage(a, 0, 0, W, H); g.globalAlpha = st.fade; g.drawImage(bb, 0, 0); g.globalAlpha = 1; }
      else if (st.mode === 'flicker') { g.drawImage(st.flick ? bb : a, 0, 0, W, H); }
      else if (st.mode === 'diff') {
        if (!diff) { diff = document.createElement('canvas'); diff.width = W; diff.height = H; const d = diff.getContext('2d'); d.drawImage(a, 0, 0, W, H); d.globalCompositeOperation = 'difference'; d.drawImage(bb, 0, 0); const id = d.getImageData(0, 0, W, H), p = id.data; for (let i = 0; i < p.length; i += 4) { const v = Math.min(255, (p[i] + p[i + 1] + p[i + 2]) * 1.6); p[i] = v; p[i + 1] = v * .35; p[i + 2] = 255 - v * .6 > 0 ? 60 : 0; p[i + 3] = 255; } d.putImageData(id, 0, 0); }
        g.drawImage(diff, 0, 0);
      }
      if (st.grid) SL.Grid.draw(g, 0, 0, W, H, { pin: st.pin });
    };
    const setPos = e => { const r = cv.getBoundingClientRect(); st.pos = clamp((e.clientX - r.left) / r.width, 0, 1); draw(); };
    let dragging = false;
    cv.addEventListener('pointerdown', e => { if (st.mode !== 'swipe') return; dragging = true; try { cv.setPointerCapture(e.pointerId); } catch (x) {} setPos(e); });
    cv.addEventListener('pointermove', e => { if (dragging) setPos(e); });
    cv.addEventListener('pointerup', () => dragging = false);
    const api = {
      set(k, v) { st[k] = v; clearInterval(st.timer); if (st.mode === 'flicker') st.timer = setInterval(() => { st.flick = !st.flick; draw(); }, 550); draw(); },
      destroy() { clearInterval(st.timer); }
    };
    api.set('mode', st.mode); return api;
  }
};

/* ----------------------------------------------------------------- scanner */
const CSS = `.sl-root{position:fixed;inset:0;z-index:1000;background:#000;color:#fff;font:15px/1.35 system-ui,-apple-system,sans-serif;display:grid;grid-template-rows:auto 1fr auto;user-select:none;-webkit-user-select:none}
.sl-view{position:relative;overflow:hidden;touch-action:none}.sl-view>video,.sl-view>canvas{position:absolute;inset:0;width:100%;height:100%;object-fit:contain}
.sl-top{display:flex;gap:6px;flex-wrap:wrap;align-items:center;padding:calc(8px + env(safe-area-inset-top)) 10px 8px;background:rgba(0,0,0,.75)}
.sl-chip{background:rgba(255,255,255,.12);border-radius:99px;padding:4px 10px;font-size:13px;white-space:nowrap}.sl-chip.ok{background:#0e7c5f}.sl-chip.warn{background:#a8431b}
.sl-x{margin-left:auto;background:none;border:0;color:#fff;font:600 15px system-ui;min-height:40px;padding:6px 10px;cursor:pointer}
.sl-guide{position:absolute;left:50%;top:14px;transform:translateX(-50%);background:rgba(0,0,0,.72);border-radius:14px;padding:8px 14px;font-weight:700;font-size:17px;text-align:center;max-width:92%;pointer-events:none}
.sl-guide.lock{background:#0e7c5f}
.sl-bot{background:rgba(0,0,0,.85);padding:10px 12px calc(12px + env(safe-area-inset-bottom));display:grid;gap:10px}
.sl-row{display:flex;gap:8px;align-items:center;flex-wrap:wrap;justify-content:center}
.sl-b{background:rgba(255,255,255,.14);color:#fff;border:0;border-radius:10px;font:600 14px system-ui;min-height:40px;padding:6px 12px;cursor:pointer}.sl-b[aria-pressed=true]{background:#fff;color:#000}
.sl-shot{width:74px;height:74px;border-radius:50%;border:0;background:none;position:relative;cursor:pointer;padding:0}
.sl-shot svg{position:absolute;inset:0}.sl-shot i{position:absolute;inset:9px;border-radius:50%;background:#fff}
.sl-meter{height:8px;border-radius:4px;background:rgba(255,255,255,.18);overflow:hidden;flex:1;min-width:80px}.sl-meter div{height:100%;background:linear-gradient(90deg,#ff4d4f,#f5b400,#19c37d);transition:width .15s}
.sl-fb{position:absolute;inset:0;display:grid;place-content:center;gap:12px;padding:24px;text-align:center;background:#0b0f14}
.sl-fb p{max-width:420px;margin:0 auto;color:#c9d3df}.sl-fb label,.sl-fb button{justify-self:center}
.sl-range{width:110px;accent-color:#fff}
@media (prefers-reduced-motion:reduce){.sl-meter div{transition:none}}`;
let cssDone = false;
const injectCSS = () => { if (cssDone) return; const s = document.createElement('style'); s.textContent = CSS; document.head.appendChild(s); cssDone = true; };

SL.Scanner = {
  /* opts: { mode:'ref'|'return', title, ref:{img, geo, pitch, roll, heading, pin, zoom}, sim:{img}, onDone(result|null) } */
  open(opts) {
    injectCSS();
    const isRet = !!(opts.mode === 'return' && opts.ref && opts.ref.img), isQR = opts.mode === 'qr';
    const S = { geo: null, ori: null,  t: null, card: null, tags: [], pin: opts.ref?.pin || null, ghost: isRet ? 'ghost' : 'off', alpha: .45, grid: !isQR, qr: null, lastQ: 0, auto: isRet, lockAt: 0, busy: false, live: false, sim: false, simT: { dx: .07, dy: -.05, rot: 5, scale: 1.1 }, closed: false, track: null, torch: false, lastA: 0, lastC: 0 };
    const stopOri = SL.Orient.start(o => { if (!o.err) S.ori = o; else S.oriErr = o.err; });
    const stopGeo = SL.Geo.watch(g => { if (!g.err) S.geo = g; else S.geoErr = g.err; });
    const root = document.createElement('div'); root.className = 'sl-root'; root.setAttribute('role', 'dialog'); root.setAttribute('aria-label', opts.title || 'Camera');
    root.innerHTML = `<div class="sl-top"><span class="sl-chip" data-c="gps">GPS…</span><span class="sl-chip" data-c="lvl">Level…</span><span class="sl-chip" data-c="card">Card: looking</span><button class="sl-x" data-a="close">Close</button></div>
    <div class="sl-view"><video playsinline muted></video><canvas class="sl-feed" hidden></canvas><canvas class="sl-ov"></canvas><div class="sl-guide">Starting camera…</div></div>
    <div class="sl-bot"><div class="sl-row"><div class="sl-meter" title="Match"><div style="width:0%"></div></div><span data-c="q" style="min-width:4em;text-align:right">${isRet ? 'Lined up 0%' : ''}</span></div>
    <div class="sl-row">${isRet ? `<button class="sl-b" data-a="ghost" aria-pressed="true">Ghost</button><button class="sl-b" data-a="edges" aria-pressed="false">Edges</button><input class="sl-range" type="range" min="10" max="90" value="45" aria-label="Ghost strength">` : ''}<button class="sl-b" data-a="grid" aria-pressed="true">Grid</button>${isRet ? '<button class="sl-b" data-a="auto" aria-pressed="true">Auto-capture</button>' : ''}<button class="sl-b" data-a="torch" hidden>Light</button></div>
    <div class="sl-row" data-r="sim" hidden><button class="sl-b" data-a="rl">Rotate left</button><button class="sl-b" data-a="zo">Smaller</button><button class="sl-b" data-a="zi">Bigger</button><button class="sl-b" data-a="rr">Rotate right</button></div>
    <div class="sl-row"><label class="sl-b" style="display:inline-flex;align-items:center">Use a photo<input type="file" accept="image/*" capture="environment" hidden data-a="file"></label><button class="sl-shot" data-a="shot" aria-label="Capture"><svg viewBox="0 0 74 74"><circle cx="37" cy="37" r="33" fill="none" stroke="rgba(255,255,255,.5)" stroke-width="5"/><circle data-c="ring" cx="37" cy="37" r="33" fill="none" stroke="#19c37d" stroke-width="5" stroke-dasharray="207.3" stroke-dashoffset="207.3" transform="rotate(-90 37 37)"/></svg><i></i></button><span style="width:94px;font-size:12px;color:#c9d3df" data-c="hint">${isQR ? 'Hold the label inside the frame' : isRet ? 'Line up the ghost' : 'Tap the grid to pin the target'}</span></div></div>`;
    document.body.appendChild(root);
    const $ = s => root.querySelector(s), video = $('video'), feed = $('.sl-feed'), ov = $('.sl-ov'), view = $('.sl-view'), guide = $('.sl-guide');
    let refImg = null, refP = null, refEdges = null, simImg = null;
    const ready = (async () => {
      if (isRet) { refImg = typeof opts.ref.img === 'string' ? await SL.img(opts.ref.img) : opts.ref.img; refP = SL.Align.prep(refImg, 128); refEdges = SL.edges(refImg); }
      if (opts.sim && opts.sim.img) simImg = typeof opts.sim.img === 'string' ? await SL.img(opts.sim.img) : opts.sim.img;
    })();
    const finish = r => { if (S.closed) return; S.closed = true; stopOri(); stopGeo(); try { (video.srcObject?.getTracks() || []).forEach(t => t.stop()); } catch (e) {} root.remove(); opts.onDone && opts.onDone(r); };
    const src = () => S.sim ? feed : video;
    const frameRect = () => { const [w, h] = S.sim ? [feed.width, feed.height] : [video.videoWidth, video.videoHeight], r = view.getBoundingClientRect(); if (!w || !h) return null; const k = Math.min(r.width / w, r.height / h); return { x: (r.width - w * k) / 2, y: (r.height - h * k) / 2, w: w * k, h: h * k, k, sw: w, sh: h }; };
    const fallback = msg => { const fb = document.createElement('div'); fb.className = 'sl-fb'; fb.innerHTML = `<b style="font-size:19px">${msg}</b><p>Live camera, GPS and tilt run when the app is installed or opened from its own web address. Here you can take or pick a photo, or try a simulated camera.</p><label class="sl-b" style="display:inline-flex;align-items:center;min-height:48px">Take or pick a photo<input type="file" accept="image/*" capture="environment" hidden data-a="file"></label>${simImg || refImg ? '<button class="sl-b" data-a="sim" style="min-height:48px">Simulate the camera</button>' : ''}`; view.appendChild(fb); };
    // camera
    ready.then(async () => {
      try {
        if (!navigator.mediaDevices?.getUserMedia) throw new Error('none');
        const stream = await navigator.mediaDevices.getUserMedia({ audio: false, video: { facingMode: { ideal: 'environment' }, width: { ideal: 1920 }, height: { ideal: 1440 } } });
        if (S.closed) { stream.getTracks().forEach(t => t.stop()); return; }
        video.srcObject = stream; await video.play(); S.live = true; S.track = stream.getVideoTracks()[0];
        const caps = S.track.getCapabilities?.() || {};
        if (caps.torch) $('[data-a=torch]').hidden = false;
        const adv = []; if (caps.focusMode?.includes('continuous')) adv.push({ focusMode: 'continuous' }); if (opts.ref?.zoom && caps.zoom) adv.push({ zoom: clamp(opts.ref.zoom, caps.zoom.min, caps.zoom.max) });
        if (adv.length) S.track.applyConstraints({ advanced: adv }).catch(() => {});
      } catch (e) { video.hidden = true; fallback(e.name === 'NotAllowedError' || e.name === 'SecurityError' ? 'Camera is blocked here' : 'No camera available'); }
      loop();
    });
    const startSim = () => { S.sim = true; video.hidden = true; feed.hidden = false; feed.width = 1200; feed.height = 900; $('[data-r=sim]').hidden = false; root.querySelector('.sl-fb')?.remove(); if (!isRet) S.simT = isQR ? { dx: .12, dy: .06, rot: 8, scale: .55 } : SL.T.id(); };
    const decodeQR = srcEl => { if (!G.jsQR) return null; const [w, h] = dims(srcEl), k = Math.min(1, 900 / Math.max(w, h)), W = Math.round(w * k), H = Math.round(h * k), c = canvasOf('qr', W, H), g = c.getContext('2d', { willReadFrequently: true }); g.drawImage(srcEl, 0, 0, W, H); const r = G.jsQR(g.getImageData(0, 0, W, H).data, W, H, { inversionAttempts: 'attemptBoth' }); if (!r || !r.data) return null; const L = r.location; return { data: r.data, poly: [L.topLeftCorner, L.topRightCorner, L.bottomRightCorner, L.bottomLeftCorner].map(p => [p.x / k, p.y / k]) }; };
    const drawSim = () => { const g = feed.getContext('2d'), im = simImg || refImg; g.fillStyle = '#1a1a1a'; g.fillRect(0, 0, 1200, 900); if (!im) return; const t = S.simT; g.save(); g.translate(600 + t.dx * 1200, 450 + t.dy * 1200); g.rotate(t.rot * D2R); g.scale(t.scale, t.scale); g.translate(-600, -450); g.drawImage(im, 0, 0, 1200, 900); g.restore(); };
    // analysis
    const analyse = now => {
      const s = src(); if (S.busy || !(S.live || S.sim)) return;
      S.busy = true;
      try {
        if (isQR && now - S.lastQ > 180 && !S.qr) { S.lastQ = now; const q = decodeQR(s); if (q) { S.qr = q; setTimeout(() => finish({ code: q.data, how: S.sim ? 'sim' : 'live', t: new Date().toISOString(), geo: S.geo }), 450); } }
        if (isRet && now - S.lastA > 140) { S.lastA = now; const p = SL.Align.prep(s, 128), r = SL.Align.register(refP, p, S.t && S.t.q > .25 ? { near: S.t, refine: false } : { refine: false });
          // low-texture scenes (grass, gravel) give weak peaks; a match that holds steady frame after frame is still a real match
          const H = S.hist = (S.hist || []).concat([r]).slice(-4), steady = H.length === 4 && H.every(h => Math.abs(h.dx - r.dx) < .008 && Math.abs(h.dy - r.dy) < .008 && Math.abs(h.rot - r.rot) < 1.1 && Math.abs(h.scale - r.scale) < .015);
          if (steady && r.peak > .2) r.q = Math.max(r.q, .45 + Math.min(.4, (r.peak - .2) * 2)); S.t = r; }
        if (!isQR && now - S.lastC > 420 && SL.Card.ready()) { S.lastC = now; const d = SL.Card.detect(s, 900); S.card = d && d.card; S.tags = d ? d.tags : [];
          if (isRet && S.card && opts.ref.card) { const [w, h] = dims(s), k = w / opts.ref.w, Hr = SL.H.mul(opts.ref.card.H, [1 / k, 0, 0, 0, 1 / k, 0, 0, 0, 1]), Hc = SL.Align.fromCards(Hr, S.card.H); if (Hc) { const tc = SL.T.fromH(Hc, w, h); S.tc = { ...tc, q: 1, at: now }; } } }
        if (S.tc && now - S.tc.at < 900) S.t = { ...S.tc, q: Math.max(S.t?.q || 0, .9) };
      } catch (e) {}
      S.busy = false;
    };
    const lvl = () => {
      if (S.sim) return { ok: true, text: 'Level: simulated' };
      if (!S.ori) return { ok: null, text: S.oriErr || 'Level…' };
      if (isRet && opts.ref.pitch != null) { const dp = S.ori.pitch - opts.ref.pitch, dr = S.ori.roll - opts.ref.roll, off = Math.hypot(dp, dr); return { ok: off < 3, off, text: off < 3 ? 'Angle matches' : `Angle off ${off.toFixed(0)}°` }; }
      return { ok: null, text: `Tilt ${S.ori.pitch.toFixed(0)}° / ${S.ori.roll.toFixed(0)}°` };
    };
    const gps = () => {
      if (S.sim) return { ok: true, text: 'GPS: simulated, at the spot' };
      if (!S.geo) return { ok: null, text: S.geoErr || 'GPS…' };
      if (isRet && opts.ref.geo) { const d = SL.Geo.dist(S.geo, opts.ref.geo), b = SL.Geo.bearing(S.geo, opts.ref.geo); return { ok: d <= Math.max(5, S.geo.acc), d, text: d < Math.max(4, S.geo.acc) ? `At the spot (±${SL.dist(S.geo.acc)})` : `${SL.dist(d)} ${SL.Geo.compass(b)} · ±${SL.dist(S.geo.acc)}` }; }
      return { ok: S.geo.acc < 10, text: `GPS ±${SL.dist(S.geo.acc)}` };
    };
    const chip = (k, o) => { const e = $(`[data-c=${k}]`); e.textContent = o.text; e.className = 'sl-chip' + (o.ok === true ? ' ok' : o.ok === false ? ' warn' : ''); };
    const render = now => {
      const fr = frameRect(), g = ov.getContext('2d'), r = view.getBoundingClientRect(), dpr = Math.min(2, devicePixelRatio || 1);
      if (ov.width !== Math.round(r.width * dpr)) { ov.width = Math.round(r.width * dpr); ov.height = Math.round(r.height * dpr); }
      g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, r.width, r.height);
      if (!fr) return;
      if (isRet && S.ghost !== 'off') { g.globalAlpha = S.ghost === 'edges' ? .95 : S.alpha; g.drawImage(S.ghost === 'edges' ? refEdges : refImg, fr.x, fr.y, fr.w, fr.h); g.globalAlpha = 1; }
      if (S.grid) SL.Grid.draw(g, fr.x, fr.y, fr.w, fr.h, { pin: S.pin });
      if (S.card) { g.strokeStyle = '#19c37d'; g.lineWidth = 3; g.beginPath(); S.card.quad.forEach((p, i) => { const x = fr.x + p[0] * fr.k, y = fr.y + p[1] * fr.k; i ? g.lineTo(x, y) : g.moveTo(x, y); }); g.closePath(); g.stroke(); }
      if (isQR) { const bw = Math.min(fr.w, fr.h) * .62, bx = fr.x + (fr.w - bw) / 2, by = fr.y + (fr.h - bw) / 2; g.strokeStyle = S.qr ? '#19c37d' : 'rgba(255,255,255,.85)'; g.lineWidth = 4; const L = bw / 5; [[bx, by, 1, 1], [bx + bw, by, -1, 1], [bx + bw, by + bw, -1, -1], [bx, by + bw, 1, -1]].forEach(([x, y, sx, sy]) => { g.beginPath(); g.moveTo(x + sx * L, y); g.lineTo(x, y); g.lineTo(x, y + sy * L); g.stroke(); }); if (S.qr) { g.fillStyle = 'rgba(25,195,125,.35)'; g.beginPath(); S.qr.poly.forEach((p, i) => { const x = fr.x + p[0] * fr.k, y = fr.y + p[1] * fr.k; i ? g.lineTo(x, y) : g.moveTo(x, y); }); g.closePath(); g.fill(); } }
      S.tags.forEach(t => { g.strokeStyle = '#f5b400'; g.lineWidth = 3; g.beginPath(); t.corners.forEach((p, i) => { const x = fr.x + p[0] * fr.k, y = fr.y + p[1] * fr.k; i ? g.lineTo(x, y) : g.moveTo(x, y); }); g.closePath(); g.stroke(); });
      // guidance
      const L = lvl(), P = gps(); chip('lvl', L); chip('gps', P);
      if (isQR) chip('card', { ok: !!S.qr, text: S.qr ? 'Code found' : 'Looking for a code' }); else chip('card', S.card ? { ok: S.card.tilt < 8, text: `${S.card.kind === 'board' ? 'Board' : 'Card'} ${String(S.card.serial).padStart(3, '0')} · ${S.card.tilt.toFixed(0)}° off square` } : S.tags.length ? { ok: true, text: `Spot tag ${S.tags[0].id}` } : { ok: null, text: 'Card: not in view' });
      let text = '', lock = false;
      if (isQR) { text = S.qr ? 'Got it: ' + (S.qr.data.split('|')[1] || S.qr.data) : (opts.qrPrompt || 'Point at a label'); lock = !!S.qr; $('.sl-meter div').style.width = S.qr ? '100%' : '0%'; }
      else if (isRet) {
        const t = S.t; const q = t ? Math.round(t.q * clamp(1 - Math.max(Math.hypot(t.dx, t.dy) / .12, Math.abs(t.rot) / 12, Math.abs(t.scale - 1) / .2), 0, 1) * 100) : 0;
        $('.sl-meter div').style.width = q + '%'; $('[data-c=q]').textContent = `Lined up ${q}%`;
        if (P.ok === false && !S.sim && P.d > 25) text = `Head ${P.text.split('·')[0].trim()}`;
        else if (!t || t.q < .12) text = 'Find the spot: line up the ghost';
        else { const adv = SL.Align.advice(t, S.sim); if (adv.length) text = adv[0].text + (adv[1] ? ' · ' + adv[1].text.toLowerCase() : ''); else if (L.ok === false) text = 'Match the phone angle: ' + L.text.toLowerCase(); else { text = 'Locked. Hold still'; lock = true; }
          // arrow toward the fix
          if (adv.length && adv[0].k === 'm') { const cx = fr.x + fr.w / 2, cy = fr.y + fr.h / 2, sx = S.sim ? -1 : 1, ax = sx * t.dx, ay = sx * t.dy, m = Math.hypot(ax, ay) || 1, L2 = Math.min(90, 40 + m * 600); g.strokeStyle = '#ffd60a'; g.fillStyle = '#ffd60a'; g.lineWidth = 6; g.beginPath(); g.moveTo(cx, cy); g.lineTo(cx + ax / m * L2, cy + ay / m * L2); g.stroke(); g.beginPath(); g.arc(cx + ax / m * L2, cy + ay / m * L2, 9, 0, 7); g.fill(); } }
      } else {
        const q = (S.card ? 50 : 0) + (L.ok !== false ? 25 : 0) + (S.pin ? 25 : 0); $('.sl-meter div').style.width = q + '%';
        const need = opts.needCard !== false;
        text = need && !S.card ? 'Put the card flat beside the target' : S.card && S.card.tilt > 8 ? `Face it square: ${S.card.tilt.toFixed(0)}° off` : S.tags.length && opts.tagMode ? `Spot tag ${S.tags[0].id} found. Capture` : !S.pin && !opts.tagMode ? 'Tap the grid on the target' : 'Ready. Capture';
        lock = (!need || S.card) && (!S.card || S.card.tilt <= 8) && (!!S.pin || !!opts.tagMode);
      }
      guide.textContent = text; guide.classList.toggle('lock', lock);
      const ring = $('[data-c=ring]');
      if (isRet && S.auto && lock) { S.lockAt = S.lockAt || now; const pr = clamp((now - S.lockAt) / 800, 0, 1); ring.setAttribute('stroke-dashoffset', 207.3 * (1 - pr)); if (pr >= 1) capture('auto'); }
      else { S.lockAt = 0; ring.setAttribute('stroke-dashoffset', 207.3); }
    };
    const loop = now => { if (S.closed) return; if (S.sim) drawSim(); analyse(now || 0); render(now || 0); requestAnimationFrame(loop); };
    const result = async (img, how, extra = {}) => {
      guide.textContent = 'Processing…';
      const out = { img: SL.toJPEG(img, 1600), how, t: new Date().toISOString(), geo: extra.geo || (S.sim ? (opts.ref?.geo || opts.sim?.geo || null) : S.geo), pitch: S.ori?.pitch ?? opts.ref?.pitch ?? null, roll: S.ori?.roll ?? opts.ref?.roll ?? null, heading: S.ori?.heading ?? null, zoom: S.track?.getSettings?.().zoom ?? null, pin: S.pin, exif: extra.exif || null };
      const im = await SL.img(out.img); out.w = im.naturalWidth; out.h = im.naturalHeight;
      const d = SL.Card.detect(im, 1400); out.card = d && d.card ? { kind: d.card.kind, serial: d.card.serial, markers: d.card.markers, H: d.card.H, tilt: d.card.tilt, quad: d.card.quad, mmPerPx: d.card.mmPerPx } : null; out.tags = d ? d.tags.map(t => t.id) : [];
      if (isRet) out.align = SL.Align.still(refImg, im, opts.ref.card ? { H: SL.H.mul(opts.ref.card.H, [opts.ref.w / out.w, 0, 0, 0, opts.ref.w / out.w, 0, 0, 0, 1]) } : null, out.card);
      finish(out);
    };
    const capture = how => { if (S.capturing) return; S.capturing = true; const s = src(); const c = document.createElement('canvas'); const [w, h] = dims(s); c.width = w; c.height = h; c.getContext('2d').drawImage(s, 0, 0); result(c, how); };
    // events
    root.addEventListener('click', e => {
      const b = e.target.closest('[data-a]'); if (!b) return; const a = b.dataset.a;
      if (a === 'close') finish(null);
      if (a === 'shot') { if (isQR) { const q = (S.live || S.sim) && decodeQR(src()); if (q) { S.qr = q; finish({ code: q.data, how: 'manual', t: new Date().toISOString(), geo: S.geo }); } else guide.textContent = 'No code in view yet'; } else if (S.live || S.sim) capture('manual'); }
      if (a === 'ghost' || a === 'edges') { S.ghost = S.ghost === a ? 'off' : a; root.querySelectorAll('[data-a=ghost],[data-a=edges]').forEach(x => x.setAttribute('aria-pressed', S.ghost === x.dataset.a)); }
      if (a === 'grid') { S.grid = !S.grid; b.setAttribute('aria-pressed', S.grid); }
      if (a === 'auto') { S.auto = !S.auto; b.setAttribute('aria-pressed', S.auto); }
      if (a === 'torch' && S.track) { S.torch = !S.torch; S.track.applyConstraints({ advanced: [{ torch: S.torch }] }).catch(() => {}); b.setAttribute('aria-pressed', S.torch); }
      if (a === 'sim') startSim();
      if (a === 'rl') S.simT.rot -= 2; if (a === 'rr') S.simT.rot += 2; if (a === 'zi') S.simT.scale *= 1.04; if (a === 'zo') S.simT.scale /= 1.04;
    });
    root.addEventListener('input', e => { if (e.target.classList.contains('sl-range')) S.alpha = e.target.value / 100; });
    root.addEventListener('change', async e => {
      if (e.target.dataset.a !== 'file' || !e.target.files[0]) return;
      const f = e.target.files[0], buf = await f.arrayBuffer().catch(() => null), exif = buf ? SL.Exif.read(buf) : {};
      if (isQR) { const url2 = URL.createObjectURL(f); try { const im = await SL.img(url2), q = decodeQR(im); if (q) finish({ code: q.data, how: 'photo', t: new Date().toISOString(), geo: S.geo }); else guide.textContent = 'No code found in that photo. Try closer.'; } catch (x) { guide.textContent = 'That photo could not be opened'; } finally { URL.revokeObjectURL(url2); } return; }
      const url = URL.createObjectURL(f); try { await ready; const im = await SL.img(url); await result(im, 'photo', { exif, geo: exif.lat != null ? { lat: exif.lat, lon: exif.lon, acc: null, src: 'photo' } : null }); } catch (x) { guide.textContent = 'That photo could not be opened'; } finally { URL.revokeObjectURL(url); }
    });
    // pinning, and simulated-camera gestures
    const pts = new Map(); let base = null;
    view.addEventListener('pointerdown', e => { if (e.target.closest('.sl-fb,button,label,input')) return; pts.set(e.pointerId, [e.clientX, e.clientY]); try { view.setPointerCapture(e.pointerId); } catch (x) {} base = { t: { ...S.simT }, pts: new Map(pts), moved: false }; });
    view.addEventListener('pointermove', e => {
      if (!pts.has(e.pointerId) || !base) return; pts.set(e.pointerId, [e.clientX, e.clientY]); const fr = frameRect(); if (!fr || !S.sim) return;
      const ids = [...pts.keys()];
      if (ids.length === 1) { const p0 = base.pts.get(ids[0]); if (!p0) return; const dx = (e.clientX - p0[0]) / fr.w, dy = (e.clientY - p0[1]) / fr.w; if (Math.hypot(dx, dy) > .005) base.moved = true; S.simT.dx = base.t.dx + dx; S.simT.dy = base.t.dy + dy; }
      else if (ids.length === 2) { const [a0, b0] = ids.map(i => base.pts.get(i)), [a1, b1] = ids.map(i => pts.get(i)); if (!a0 || !b0) return; base.moved = true; S.simT.scale = base.t.scale * Math.hypot(b1[0] - a1[0], b1[1] - a1[1]) / (Math.hypot(b0[0] - a0[0], b0[1] - a0[1]) || 1); S.simT.rot = base.t.rot + (Math.atan2(b1[1] - a1[1], b1[0] - a1[0]) - Math.atan2(b0[1] - a0[1], b0[0] - a0[0])) * R2D; }
    });
    const up = e => { if (!pts.has(e.pointerId)) return; const was = base; pts.delete(e.pointerId); if (pts.size === 0) { if (was && !was.moved && !isRet) { const fr = frameRect(), r = view.getBoundingClientRect(); if (fr) { const fx = (e.clientX - r.left - fr.x) / fr.w, fy = (e.clientY - r.top - fr.y) / fr.h; if (fx >= 0 && fx <= 1 && fy >= 0 && fy <= 1) { S.pin = { x: fx, y: fy, cell: SL.Grid.label(fx, fy) }; $('[data-c=hint]').textContent = 'Target pinned at ' + S.pin.cell; } } } base = null; } else base = { t: { ...S.simT }, pts: new Map(pts), moved: true }; };
    view.addEventListener('pointerup', up); view.addEventListener('pointercancel', up);
    view.addEventListener('wheel', e => { if (!S.sim) return; e.preventDefault(); S.simT.scale *= e.deltaY < 0 ? 1.03 : 1 / 1.03; }, { passive: false });
    const esc = e => { if (e.key === 'Escape') { removeEventListener('keydown', esc); finish(null); } }; addEventListener('keydown', esc);
    return { close: () => finish(null), _state: S, _startSim: startSim };
  }
};

G.SiteLock = SL;
})(typeof window !== 'undefined' ? window : globalThis);
