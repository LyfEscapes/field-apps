// Cylinder Chain offline cache · 2026.11.23-1
const CACHE='cylinderchain-2026.11.23-1';
const SHELL=["./", "icons/icon-192.png", "icons/icon-512.png", "index.html", "manifest.webmanifest", "privacy.html", "sitelock.js", "vendor/aruco.js", "vendor/cv.js", "vendor/jsQR.js", "vendor/qrcode.js"];
self.addEventListener('install',e=>{e.waitUntil(caches.open(CACHE).then(c=>c.addAll(SHELL)).then(()=>self.skipWaiting()))});
self.addEventListener('activate',e=>{e.waitUntil(caches.keys().then(ks=>Promise.all(ks.filter(k=>k.startsWith('cylinderchain-')&&k!==CACHE).map(k=>caches.delete(k)))).then(()=>self.clients.claim()))});
self.addEventListener('fetch',e=>{const u=new URL(e.request.url);if(e.request.method!=='GET')return;
  if(u.origin===location.origin){e.respondWith(caches.match(e.request,{ignoreSearch:true}).then(hit=>{const net=fetch(e.request).then(r=>{if(r.ok)caches.open(CACHE).then(c=>c.put(e.request,r.clone()));return r}).catch(()=>hit);return hit||net}));return}
  if(/fonts\.(googleapis|gstatic)\.com/.test(u.host)){e.respondWith(caches.open(CACHE+'-fonts').then(c=>c.match(e.request).then(hit=>hit||fetch(e.request).then(r=>{c.put(e.request,r.clone());return r}))))}
});
