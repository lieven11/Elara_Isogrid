#!/usr/bin/env python3
"""
Version 2 dashboard: adds lattice selection (isogrid 60°, triangular with variable angle, honeycomb),
optional optimization over material and lattice, a small pattern preview canvas, and finer black gridlines.
"""
from pathlib import Path
import argparse
import json

HTML_TEMPLATE = """<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>Isogrid/Honeycomb Dashboard v2</title>
  <script src=\"https://cdn.plot.ly/plotly-2.27.0.min.js\"></script>
  <style>
    body { font-family: system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin:0; padding:0; }
    .wrap { display:grid; grid-template-columns: 360px 1fr; gap:12px; height:100vh; }
    .panel { overflow:auto; padding: 12px 14px; border-right:1px solid #e5e5e5; }
    .panel h2 { margin:0 0 6px 0; font-size:18px; }
    .subhead { font-weight:600; font-size:12px; color:#222; margin:12px 0 6px; border-bottom:1px solid #eee; padding-bottom:3px; }
    .row { display:grid; grid-template-columns: 1fr 1fr; gap:8px; }
    .row3 { display:grid; grid-template-columns: 1fr 1fr 1fr; gap:8px; }
    input[type=range]{ width:100%%; }
    .muted{ color:#666; font-size:12px; }
    .badge { background:#eef; color:#224; padding:2px 6px; border-radius:3px; font-size:11px; }
    #previewBox{ display:none; margin-top:6px; }
    #previewCanvas{ width:320px; height:240px; border:1px solid #ddd; background:#fff; }
  </style>
</head>
<body>
  <div class=\"wrap\">
    <div class=\"panel\">
      <h2>Dashboard v2 <span class=\"badge\">isogrid/honeycomb</span></h2>
      <div class=\"subhead\">Achsen</div>
      <div class=\"ctrl\">
        <label>Axes (3 design parameters)</label>
        <div class=\"row3\">
          <select id=\"xAxis\"></select>
          <select id=\"yAxis\"></select>
          <select id=\"zAxis\"></select>
        </div>
        <div class=\"muted\">Choose distinct parameters for x, y, z.</div>
      </div>

      <div class=\"subhead\">Metric</div>
      <div class=\"ctrl\"><select id=\"metric\"></select></div>

      <div class=\"subhead\">Material & Lattice</div>
      <div class=\"ctrl\">
        <div class=\"row\">
          <div><label>Material</label><select id=\"material\"></select></div>
          <div><label>Lattice</label><select id=\"lattice\"></select></div>
        </div>
        <div class=\"row\"><div><label>Preview</label><input id=\"showPreview\" type=\"checkbox\" /> <span class=\"muted\">show pattern</span></div></div>
        <div id=\"previewBox\"><canvas id=\"previewCanvas\" width=\"320\" height=\"240\"></canvas></div>
      </div>

      <div class=\"subhead\">Sampling</div>
      <div class=\"ctrl\">
        <div class=\"row3\">
          <input id=\"nx\" type=\"range\" min=\"5\" max=\"50\" step=\"1\" value=\"%(nx)d\" />
          <input id=\"ny\" type=\"range\" min=\"5\" max=\"50\" step=\"1\" value=\"%(ny)d\" />
          <input id=\"nz\" type=\"range\" min=\"5\" max=\"50\" step=\"1\" value=\"%(nz)d\" />
        </div>
        <div class=\"row3 muted\"><div>x:<span id=\"nxVal\"></span></div><div>y:<span id=\"nyVal\"></span></div><div>z:<span id=\"nzVal\"></span></div></div>
      </div>

      <div class=\"subhead\">Optimizer</div>
      <div class=\"ctrl\">
        <div class=\"row\">
          <div><label>Ziel</label><select id=\"optObjective\"></select></div>
          <div><label class=\"muted\">SF_min‑Grenze</label><input id=\"sfThresh\" type=\"number\" step=\"0.05\" value=\"1.20\" style=\"width:100%%;\" /></div>
        </div>
        <div class=\"row\">
          <div><label><input type=\"checkbox\" id=\"optIncludeMat\" /> Material mitoptimieren</label></div>
          <div><label><input type=\"checkbox\" id=\"optIncludeLat\" /> Lattice mitoptimieren</label></div>
        </div>
        <div><button id=\"runOpt\">Optimize</button></div>
        <div id=\"optOut\" class=\"muted\"></div>
      </div>

      <div class=\"subhead\">Parameter‑Ranges</div>
      <div class=\"ctrl\">
        <div class=\"ctrl\"><div>b [mm]</div>
          <input id=\"bMin\" type=\"range\" min=\"1\" max=\"100\" step=\"0.5\" value=\"%(b_min_mm).1f\" />
          <input id=\"bMax\" type=\"range\" min=\"1\" max=\"100\" step=\"0.5\" value=\"%(b_max_mm).1f\" />
          <div class=\"muted\">min=<span id=\"bMinVal\"></span> max=<span id=\"bMaxVal\"></span></div>
        </div>
        <div class=\"ctrl\"><div>t [mm]</div>
          <input id=\"tMin\" type=\"range\" min=\"0.2\" max=\"20\" step=\"0.1\" value=\"%(t_min_mm).1f\" />
          <input id=\"tMax\" type=\"range\" min=\"0.2\" max=\"20\" step=\"0.1\" value=\"%(t_max_mm).1f\" />
          <div class=\"muted\">min=<span id=\"tMinVal\"></span> max=<span id=\"tMaxVal\"></span></div>
        </div>
        <div class=\"ctrl\"><div>n_theta [-]</div>
          <input id=\"nMin\" type=\"range\" min=\"10\" max=\"400\" step=\"1\" value=\"%(n_min)d\" />
          <input id=\"nMax\" type=\"range\" min=\"10\" max=\"400\" step=\"1\" value=\"%(n_max)d\" />
          <div class=\"muted\">min=<span id=\"nMinVal\"></span> max=<span id=\"nMaxVal\"></span></div>
        </div>
        <div class=\"ctrl\"><div>a [mm]</div>
          <input id=\"aMin\" type=\"range\" min=\"5\" max=\"200\" step=\"1\" value=\"%(a_min_mm).0f\" />
          <input id=\"aMax\" type=\"range\" min=\"5\" max=\"200\" step=\"1\" value=\"%(a_max_mm).0f\" />
          <div class=\"muted\">min=<span id=\"aMinVal\"></span> max=<span id=\"aMaxVal\"></span></div>
        </div>
        <div class=\"ctrl\"><div>L [m]</div>
          <input id=\"lMin\" type=\"range\" min=\"0.1\" max=\"5\" step=\"0.05\" value=\"%(L_min).2f\" />
          <input id=\"lMax\" type=\"range\" min=\"0.1\" max=\"5\" step=\"0.05\" value=\"%(L_max).2f\" />
          <div class=\"muted\">min=<span id=\"lMinVal\"></span> max=<span id=\"lMaxVal\"></span></div>
        </div>
      </div>

      <div class=\"subhead\">Fixwerte/Lasten</div>
      <div class=\"ctrl\">
        <div class=\"row\"><div><div>R = <span id=\"RVal\"></span> mm</div><input id=\"Rmm\" type=\"range\" min=\"50\" max=\"5000\" step=\"10\" value=\"%(R_mm).0f\" /></div>
          <div><div>K (Euler) = <span id=\"KVal\"></span></div><input id=\"K\" type=\"range\" min=\"0.5\" max=\"2.0\" step=\"0.05\" value=\"%(K).2f\" /></div></div>
        <div class=\"row\"><div><div>KDF [-] = <span id=\"KDFVal\"></span></div><input id=\"KDF\" type=\"range\" min=\"0.40\" max=\"1.00\" step=\"0.01\" value=\"%(KDF).2f\" /></div></div>
        <div class=\"row\"><div><div>N_req [N] = <span id=\"NVal\"></span></div><input id=\"Nreq\" type=\"range\" min=\"0\" max=\"500000\" step=\"1000\" value=\"%(N_req).0f\" /></div>
          <div><div>M_req [N·m] = <span id=\"MVal\"></span></div><input id=\"Mreq\" type=\"range\" min=\"0\" max=\"5000\" step=\"10\" value=\"%(M_req).0f\" /></div></div>
        <div class=\"row\"><div><div>T_req [N·m] = <span id=\"TVal\"></span></div><input id=\"Treq\" type=\"range\" min=\"0\" max=\"5000\" step=\"10\" value=\"%(T_req).0f\" /></div></div>
      </div>
      <div class=\"muted\">Torsion: k_torsion=0.2; KDF global als Knock‑down. Lattice‑Faktoren stark vereinfacht (Triangular: F≈2·sinθ; Honeycomb: F≈2/√3).</div>
    </div>
    <div id=\"plot\" style=\"height:100%%;\"></div>
  </div>

<script>
// Materials
const MAT_DB = {
  "AISI 304":   {E: 193e9, nu: 0.29, rho: 8000, sigma_y: 215e6},
  "AISI 316L":  {E: 193e9, nu: 0.30, rho: 8000, sigma_y: 170e6},
  "Al 6061-T6": {E: 69e9,  nu: 0.33, rho: 2700, sigma_y: 276e6},
  "Al 2024-T3": {E: 73e9,  nu: 0.33, rho: 2780, sigma_y: 325e6},
  "Al 7075-T6": {E: 72e9,  nu: 0.33, rho: 2810, sigma_y: 505e6}
};
let mat = MAT_DB["AISI 304"]; // current material
const FOS = 1.25, k_crippling = 4.0, k_torsion = 0.2;

// Lattices
const LATTICES = {
  tri_0:      {label:'Triangle 0°', kind:'tri', orient:0},
  tri_90:     {label:'Triangle 90°', kind:'tri', orient:90},
  tri_hex:    {label:'Tri-hexagonal', kind:'tri', orient:0},
  honeycomb:  {label:'Hexagonal (honeycomb)', kind:'hex'},
  square_0:   {label:'Squares 0°', kind:'sq', orient:0},
  square_45:  {label:'Squares 45°', kind:'sq', orient:45},
};

const defaultState = %(state_json)s;
let GEOM = {S:Math.sqrt(3), dz_fac:Math.sqrt(3)/2, theta:Math.PI/3, type:'tri_0', orient:0, fams:[0,Math.PI/3,-Math.PI/3], axial_scale:1.0}; // updated from state

function mm_to_m(x){ return x/1000.0 }
function geomFromState(state){
  const type = state.lattice || 'tri_0';
  let S, dz_fac, theta = Math.PI/3, orient=0, fams=[0,Math.PI/3,-Math.PI/3];
  const lat = LATTICES[type] || LATTICES.tri_0;
  orient = (lat.orient||0) * Math.PI/180.0;
  if(lat.kind==='hex'){
    // hex edges: 0°, ±60° families
    S = 2/Math.sqrt(3); dz_fac = Math.sqrt(3)/2; theta = Math.PI/3; fams=[0,Math.PI/3,-Math.PI/3];
  } else if(lat.kind==='sq') {
    // square grid: two orthogonal families
    S = 2.0; dz_fac = (type==='square_45') ? Math.SQRT1_2 : 1.0; theta = Math.PI/4; fams=[0,Math.PI/2];
  } else { // triangular families
    S = Math.sqrt(3); dz_fac = Math.sqrt(3)/2; theta = Math.PI/3; fams=[0,Math.PI/3,-Math.PI/3];
  }
  // axial scaling for EA/EI from orientation content; normalize to tri_0 sum(1 + 0.25 + 0.25) = 1.5
  const sum_c2 = fams.reduce((acc,a)=> acc + Math.cos(a+orient)**2, 0);
  const axial_scale = sum_c2 / 1.5; // tri_0 baseline
  return {S, dz_fac, theta, type, orient, fams, axial_scale};
}
function updateGeom(state){ GEOM = geomFromState(state); }

// Mechanics using GEOM
function areal_mass(b,t,a){ return mat.rho*(b*t)*(GEOM.S/Math.max(a,1e-12)); }
function total_mass(R,L,b,t,a){ return areal_mass(b,t,a)*(2*Math.PI*R*L); }
function EI_eq(R,b,t,a){ return GEOM.axial_scale * Math.PI*mat.E*(b*t)*(GEOM.S/Math.max(a,1e-12))*R**3; }
function Ncr_global(R,L,b,t,a,K){ return (Math.PI**2)*EI_eq(R,b,t,a)/(Math.max(K*L,1e-12))**2; }
function EA_ring(R,b,t,a){ return GEOM.axial_scale * Math.PI*mat.E*(b*t)*(GEOM.S/Math.max(a,1e-12))*R; }
function delta_z(a){ return GEOM.dz_fac*a; }
function Imin_rect(b,t){ return b*t**3/12; }
function euler_sigma_cr(len,r_g,K){ const lam=(K*len)/Math.max(r_g,1e-12); return (Math.PI**2)*mat.E/(lam**2); }
function crippling_sigma(b,t){ const denom=12*(1-mat.nu**2); return k_crippling*(Math.PI**2)*mat.E/denom*(t/Math.max(b,1e-12))**2; }
function local_sf(R,L,b,t,a,N_req,K){
  if(!N_req||N_req<=0) return Infinity;
  const EA=EA_ring(R,b,t,a); const eps=N_req/Math.max(EA,1e-12);
  const Imin=Imin_rect(b,t); const A=b*t; const rg=Math.sqrt(Math.max(Imin/Math.max(A,1e-12),1e-18));
  const sigma_allow = mat.sigma_y/FOS; const sigma_crip = crippling_sigma(b,t);
  let sfmin = Infinity;
  // families with angles GEOM.fams (rad) + orient
  for(const baseAng of GEOM.fams){
    const ang = baseAng + (GEOM.orient||0);
    const cos2 = Math.cos(ang)**2;
    const dem = mat.E*eps*cos2;
    // choose effective length: near-axial uses dz, others use a (approx.)
    const nearAxial = Math.abs(Math.cos(ang)) > 0.95;
    const len = nearAxial ? delta_z(a) : a;
    const se = euler_sigma_cr(len, rg, K);
    const scrit = Math.min(se, sigma_crip, sigma_allow);
    const sf = dem<=0 ? Infinity : scrit/dem;
    if(sf < sfmin) sfmin = sf;
  }
  return sfmin;
}
function GJ_eq(R,b,t,a){ const G=mat.E/(2*(1+mat.nu)); return k_torsion*Math.PI*G*(b*t)*(GEOM.S/Math.max(a,1e-12))*R**3; }
function sf_bending(R,b,t,a,M_req){ if(!M_req||M_req<=0) return NaN; const sigma_allow=mat.sigma_y/FOS; const M_cap=sigma_allow*EI_eq(R,b,t,a)/Math.max(R,1e-12); return M_cap/Math.max(M_req,1e-12); }
function sf_torsion(R,b,t,a,T_req){ if(!T_req||T_req<=0) return NaN; const tau_allow=0.58*(mat.sigma_y/FOS); const J=GJ_eq(R,b,t,a); const tau=T_req*Math.max(R,0)/Math.max(J,1e-18); return tau_allow/Math.max(tau,1e-18); }

const METRICS = {
  ncr_global: {label:'Ncr [N]', fn:(R,L,b,t,a,K,KDF,N,M,T)=> KDF*Ncr_global(R,L,b,t,a,K)},
  mass: {label:'Masse [kg]', fn:(R,L,b,t,a,K,KDF)=> total_mass(R,L,b,t,a)},
  areal_mass: {label:'Arealmasse [kg/m²]', fn:(R,L,b,t,a)=> areal_mass(b,t,a)},
  ncr_over_m: {label:'Ncr/m [m/s²]', fn:(R,L,b,t,a,K,KDF)=> (KDF*Ncr_global(R,L,b,t,a,K))/Math.max(total_mass(R,L,b,t,a),1e-12)},
  sf_min: {label:'SF_min [-]', fn:(R,L,b,t,a,K,KDF,N)=> Math.min(((KDF*Ncr_global(R,L,b,t,a,K))/Math.max(N||1e-12,1e-12)), local_sf(R,L,b,t,a,N,K))},
  sf_bending: {label:'SF Biegung [-]', fn:(R,L,b,t,a,K,KDF,N,M)=> sf_bending(R,b,t,a,M)},
  sf_torsion: {label:'SF Torsion [-]', fn:(R,L,b,t,a,K,KDF,N,M,T)=> sf_torsion(R,b,t,a,T)}
};

const AX = ['b','t','n_theta','a','L'];
const AX_LABEL = {b:'b [mm]', t:'t [mm]', n_theta:'n_theta [-]', a:'a [mm]', L:'L [m]'};

function linspace(minv,maxv,n){ if(n<=1) return [minv]; const step=(maxv-minv==0)?1:(maxv-minv)/(n-1); return Array.from({length:n},(_,i)=>minv+i*step); }

function n_theta_from_a(R,a){ return Math.max(1, Math.round(2*Math.PI*R/Math.max(a,1e-12))); }
function a_from_n_theta(R,n){ return 2*Math.PI*R/Math.max(1,n); }

function buildSurface(state){
  updateGeom(state);
  const R=mm_to_m(state.R_mm), K=state.K, KDF=state.KDF, N=state.N_req, M=state.M_req, T=state.T_req;
  const ax=state.xAxis, ay=state.yAxis; const nx=state.nx, ny=state.ny;
  function rng(name, n){
    if(name==='b') return linspace(state.b_min_mm, state.b_max_mm, n).map(mm_to_m);
    if(name==='t') return linspace(state.t_min_mm, state.t_max_mm, n).map(mm_to_m);
    if(name==='n_theta') return linspace(state.n_min, state.n_max, n).map(v=>Math.round(v));
    if(name==='a') return linspace(state.a_min_mm, state.a_max_mm, n).map(mm_to_m);
    if(name==='L') return linspace(state.L_min, state.L_max, n);
    return [];
  }
  const Xs=rng(ax,nx), Ys=rng(ay,ny); const Z=[];
  for(let j=0;j<Ys.length;j++){
    const row=[];
    for(let i=0;i<Xs.length;i++){
      let b=mm_to_m(state.b_fix_mm), t=mm_to_m(state.t_fix_mm), nth=Math.round(state.n_fix), a=mm_to_m(state.a_fix_mm), L=state.L_fix;
      function apply(nm,val){ if(nm==='b') b=val; if(nm==='t') t=val; if(nm==='n_theta') nth=val; if(nm==='a') a=val; if(nm==='L') L=val; }
      apply(ax, Xs[i]); apply(ay, Ys[j]);
      if(state.use_a_priority){ nth=n_theta_from_a(R,a);} else { a=a_from_n_theta(R,nth);}    
      const val=METRICS[state.metric].fn(R,L,b,t,a,K,KDF,N,M,T);
      row.push(Number.isFinite(val)?val:NaN);
    }
    Z.push(row);
  }
  const Xd = Xs.map(v=> (ax==='b'||ax==='t'||ax==='a')? v*1000.0 : v);
  const Yd = Ys.map(v=> (ay==='b'||ay==='t'||ay==='a')? v*1000.0 : v);
  return {X:Xd, Y:Yd, Z};
}

function axisLabel(name,b,t,nth,a,L){ if(name==='b'||name==='t'||name==='a') return (name==='a'?a*1000.0: (name==='b'?b*1000.0:t*1000.0)); if(name==='n_theta') return nth; return L; }

function drawPreviewIfNeeded(state){
  const cb = document.getElementById('showPreview');
  const box = document.getElementById('previewBox');
  if(!cb || !box) return;
  if(!cb.checked){ box.style.display='none'; return; }
  box.style.display='block';
  const c = document.getElementById('previewCanvas');
  const ctx = c.getContext('2d');
  ctx.clearRect(0,0,c.width,c.height);
  // Determine current parameters (use slice for a/n_theta if needed)
  const R=mm_to_m(state.R_mm);
  let a = mm_to_m(state.a_fix_mm), nth = Math.round(state.n_fix);
  const az = state.zAxis; // use slice
  const zdisp = state.sliceDisp;
  function fromDisplay(name, val){ if(name==='b'||name==='t'||name==='a') return mm_to_m(val); if(name==='n_theta') return Math.round(val); return val; }
  if(az==='a') a = fromDisplay('a', zdisp);
  if(az==='n_theta') nth = fromDisplay('n_theta', zdisp);
  if(state.use_a_priority){ nth = n_theta_from_a(R,a);} else { a = a_from_n_theta(R,nth); }
  const a_mm = a*1000.0;
  const b_mm = state.b_fix_mm; const t_mm = state.t_fix_mm;
  // scale
  const margin=10; const cell= Math.max(10, Math.min((c.width-2*margin)/3, (c.height-2*margin)/3));
  const s = cell / a_mm; // px per mm
  ctx.lineWidth = Math.max(1, t_mm*s);
  ctx.strokeStyle = '#333';
  // draw pattern approx: three families of lines for tri; hex grid for honeycomb
  const type = state.lattice; const orient = GEOM.orient || 0; const theta = GEOM.theta;
  if(type==='honeycomb'){
    // draw hexagons with side a_mm (tiling)
    const s_px = a_mm*s;
    const h = Math.sin(Math.PI/3)*s_px; // equilateral height
    const rowH = 2*h; const colW = 1.5*s_px;
    for(let row=0; row<5; row++){
      for(let col=0; col<5; col++){
        const cx = margin + col*colW + ((row%%2)? 0.75*s_px : 0);
        const cy = margin + row*rowH + h;
        const pts = [
          [cx - 0.5*s_px, cy - h], [cx + 0.5*s_px, cy - h], [cx + s_px, cy],
          [cx + 0.5*s_px, cy + h], [cx - 0.5*s_px, cy + h], [cx - s_px, cy]
        ];
        ctx.beginPath(); ctx.moveTo(pts[0][0], pts[0][1]); for(let i=1;i<pts.length;i++){ ctx.lineTo(pts[i][0], pts[i][1]); } ctx.closePath(); ctx.stroke();
      }
    }
  } else if(type==='square_0' || type==='square_45'){
    const s_px = a_mm*s;
    const ang1 = (type==='square_45') ? (Math.PI/4) : 0;
    const ang2 = ang1 + Math.PI/2;
    const angles = [ang1+orient, ang2+orient];
    const centerX = c.width/2, centerY = c.height/2;
    const L = Math.hypot(c.width, c.height);
    angles.forEach(phi=>{
      const ux = Math.cos(phi), uy = Math.sin(phi);
      for(let k=-10;k<20;k++){
        const px = centerX + k*s_px*(-uy); // shift perpendicular by cell size
        const py = centerY + k*s_px*(ux);
        ctx.beginPath(); ctx.moveTo(px - ux*L, py - uy*L); ctx.lineTo(px + ux*L, py + uy*L); ctx.stroke();
      }
    });
  } else if(type==='tri_hex'){
    // Kagome-like: segmented lines to suggest triangles + hexes
    const s_px = a_mm*s; const angs = [0, +theta, -theta].map(a=>a+orient);
    const centerX = c.width/2, centerY = c.height/2; const L = Math.hypot(c.width, c.height);
    ctx.lineWidth = Math.max(1, t_mm*s);
    angs.forEach(phi=>{
      const ux=Math.cos(phi), uy=Math.sin(phi); const px0=centerX, py0=centerY;
      const step = s_px; const seg = s_px*0.8; // 20%% gap
      for(let k=-10;k<11;k++){
        const px = px0 + k*step*(-uy); const py = py0 + k*step*(ux);
        const N = Math.ceil(L/seg);
        for(let n=-N;n<=N;n++){
          const t0 = n*(seg+step*0.2) - seg/2, t1 = t0 + seg;
          ctx.beginPath(); ctx.moveTo(px + ux*t0, py + uy*t0); ctx.lineTo(px + ux*t1, py + uy*t1); ctx.stroke();
        }
      }
    });
  } else { // triangular families
    const s_px = a_mm*s; const angs = [0, +theta, -theta].map(a=>a+orient);
    const centerX = c.width/2, centerY = c.height/2; const L = Math.hypot(c.width, c.height);
    angs.forEach(phi=>{
      const ux=Math.cos(phi), uy=Math.sin(phi); const px0=centerX, py0=centerY;
      for(let k=-12;k<13;k++){
        const px = px0 + k*s_px*(-uy); const py= py0 + k*s_px*(ux);
        ctx.beginPath(); ctx.moveTo(px - ux*L, py - uy*L); ctx.lineTo(px + ux*L, py + uy*L); ctx.stroke();
      }
    });
  }
}

function mount(){
  const state = Object.assign({}, defaultState);
  updateGeom(state);
  // populate selects
  const metricSel = document.getElementById('metric'); Object.keys(METRICS).forEach(k=>{ const o=document.createElement('option'); o.value=k; o.text=METRICS[k].label; metricSel.add(o); }); metricSel.value=state.metric;
  ['xAxis','yAxis','zAxis'].forEach((id,i)=>{ const sel=document.getElementById(id); AX.forEach(ax=>{const o=document.createElement('option'); o.value=ax; o.text=AX_LABEL[ax]; sel.add(o);}); sel.value=state[["xAxis","yAxis","zAxis"][i]]; });
  const matSel=document.getElementById('material'); Object.keys(MAT_DB).forEach(k=>{const o=document.createElement('option'); o.value=k; o.text=k; matSel.add(o);}); matSel.value=state.material; mat=MAT_DB[state.material];
  const latSel=document.getElementById('lattice'); Object.keys(LATTICES).forEach(k=>{ const o=document.createElement('option'); o.value=k; o.text=LATTICES[k].label; latSel.add(o);}); latSel.value=state.lattice;
  matSel.addEventListener('change', (e)=>{ state.material=e.target.value; mat=MAT_DB[state.material]; rerender(); });
  latSel.addEventListener('change', (e)=>{ state.lattice=e.target.value; updateGeom(state); rerender(); drawPreviewIfNeeded(state); });
  document.getElementById('showPreview').addEventListener('change', ()=> drawPreviewIfNeeded(state));
  document.getElementById('metric').addEventListener('change',(e)=>{ state.metric=e.target.value; rerender(); });
  ['xAxis','yAxis','zAxis'].forEach(id=> document.getElementById(id).addEventListener('change',(e)=>{ state[id]=e.target.value; rerender(); }));
  function bind(id, spanId, prop){ const el=document.getElementById(id); const span=document.getElementById(spanId); const sync=()=>{ span.textContent=(+el.value).toString(); state[prop]=+el.value; rerender(); drawPreviewIfNeeded(state); }; el.addEventListener('input', sync); span.textContent=el.value; state[prop]=+el.value; }
  bind('nx','nxVal','nx'); bind('ny','nyVal','ny'); bind('nz','nzVal','nz');
  bind('bMin','bMinVal','b_min_mm'); bind('bMax','bMaxVal','b_max_mm'); bind('tMin','tMinVal','t_min_mm'); bind('tMax','tMaxVal','t_max_mm');
  bind('nMin','nMinVal','n_min'); bind('nMax','nMaxVal','n_max'); bind('aMin','aMinVal','a_min_mm'); bind('aMax','aMaxVal','a_max_mm');
  bind('lMin','lMinVal','L_min'); bind('lMax','lMaxVal','L_max'); bind('Rmm','RVal','R_mm'); bind('K','KVal','K'); bind('KDF','KDFVal','KDF');
  bind('Nreq','NVal','N_req'); bind('Mreq','MVal','M_req'); bind('Treq','TVal','T_req');

  // slice controls
  const sliceAxisSel = document.createElement('select'); ['b','t','n_theta','a','L'].forEach(ax=>{ const o=document.createElement('option'); o.value=ax; o.text=AX_LABEL[ax]; sliceAxisSel.add(o); });
  // embed into DOM? v1 had hidden slice UI; in v2 we implicitly use zAxis as slice
  state.use_a_priority = true;

  const plotDiv=document.getElementById('plot');
  // init slice value to mid of z-axis range
  function axisRangeDisplay(st, name){
    if(name==='b') return [st.b_min_mm, st.b_max_mm];
    if(name==='t') return [st.t_min_mm, st.t_max_mm];
    if(name==='a') return [st.a_min_mm, st.a_max_mm];
    if(name==='n_theta') return [st.n_min, st.n_max];
    if(name==='L') return [st.L_min, st.L_max];
    return [0,1];
  }
  { const [zmin,zmax]=axisRangeDisplay(state, state.zAxis); state.sliceDisp = (zmin+zmax)/2; }
  let rerenderTimer=null;
  function rerender(){ if(rerenderTimer) clearTimeout(rerenderTimer); rerenderTimer = setTimeout(_rerender, 40); }
  function _rerender(){
    updateGeom(state);
    const grid=buildSurface(state);
    const ctitle = METRICS[state.metric].label;
    const layout = { title: `Surface — ${ctitle}`, scene:{ xaxis:{title:AX_LABEL[state.xAxis]}, yaxis:{title:AX_LABEL[state.yAxis]}, zaxis:{title:ctitle} }, margin:{l:0,r:0,b:0,t:38}, height:window.innerHeight, uirevision:'v2' };
    if (plotDiv && plotDiv._fullLayout && plotDiv._fullLayout.scene && plotDiv._fullLayout.scene.camera) { layout.scene.camera = plotDiv._fullLayout.scene.camera; }
    // compute fine contour steps (1/5 of sampling step)
    const dx = (grid.X.length>1) ? (grid.X[1]-grid.X[0]) : 1;
    const dy = (grid.Y.length>1) ? (grid.Y[1]-grid.Y[0]) : 1;
    const data=[{ type:'surface', x:grid.X, y:grid.Y, z:grid.Z, colorscale:'Viridis', showscale:true,
      contours:{ x:{show:true, color:'black', width:1, start:grid.X[0], end:grid.X[grid.X.length-1], size: dx/5},
                 y:{show:true, color:'black', width:1, start:grid.Y[0], end:grid.Y[grid.Y.length-1], size: dy/5}, z:{show:false} } }];
    if (state.bestMarker){ data.push({type:'scatter3d', mode:'markers', x:[state.bestMarker.x], y:[state.bestMarker.y], z:[state.bestMarker.z], marker:{size:6, color:'red'}, name:'Best'}); }
    Plotly.react(plotDiv, data, layout, {responsive:true});
    drawPreviewIfNeeded(state);
  }
  _rerender();

  // Optimizer
  const objectives = { min_mass:'Masse minimieren (SF_min ≥ Ziel)', max_ncr_over_m:'Ncr/m maximieren (SF_min ≥ Ziel)' };
  const objSel = document.getElementById('optObjective'); Object.keys(objectives).forEach(k=>{const o=document.createElement('option'); o.value=k; o.text=objectives[k]; objSel.add(o);}); objSel.value='min_mass';
  document.getElementById('runOpt').addEventListener('click', ()=> runOptimize());

  function runOptimize(){
    const sfReq = parseFloat(document.getElementById('sfThresh').value || '1.2');
    const obj = document.getElementById('optObjective').value;
    const includeMat = !!document.getElementById('optIncludeMat').checked;
    const includeLat = !!document.getElementById('optIncludeLat').checked;
    const mats = includeMat ? Object.keys(MAT_DB) : [state.material];
    const lats = includeLat ? Object.keys(LATTICES) : [state.lattice];
    function rng(name, n){ if(name==='b') return linspace(state.b_min_mm, state.b_max_mm, state.nx).map(mm_to_m);
      if(name==='t') return linspace(state.t_min_mm, state.t_max_mm, state.ny).map(mm_to_m);
      if(name==='n_theta') return linspace(state.n_min, state.n_max, state.nz).map(v=>Math.round(v));
      if(name==='a') return linspace(state.a_min_mm, state.a_max_mm, state.nz).map(mm_to_m);
      if(name==='L') return linspace(state.L_min, state.L_max, state.nz); return []; }
    const ax=state.xAxis, ay=state.yAxis, az=state.zAxis;
    const Xs=rng(ax,state.nx), Ys=rng(ay,state.ny), Zs=rng(az,state.nz);
    let best=null, bestObj=Infinity;
    for(const mname of mats){
      const matSaved = mat; const stateMatSaved = state.material; mat = MAT_DB[mname]; state.material = mname;
      for(const lat of lats){
        const latSaved=state.lattice; state.lattice=lat; updateGeom(state);
          for(let k=0;k<Zs.length;k++){
            for(let j=0;j<Ys.length;j++){
              for(let i=0;i<Xs.length;i++){
                const R=mm_to_m(state.R_mm), K=state.K, KDF=state.KDF, N=state.N_req, M=state.M_req, T=state.T_req;
                let b=mm_to_m(state.b_fix_mm), t=mm_to_m(state.t_fix_mm), nth=Math.round(state.n_fix), a=mm_to_m(state.a_fix_mm), L=state.L_fix;
                function put(name,val){ if(name==='b') b=val; if(name==='t') t=val; if(name==='n_theta') nth=val; if(name==='a') a=val; if(name==='L') L=val; }
                put(ax,Xs[i]); put(ay,Ys[j]); put(az,Zs[k]);
                if(state.use_a_priority){ nth = n_theta_from_a(R,a);} else { a = a_from_n_theta(R,nth);}                
                const ncr = KDF * Ncr_global(R,L,b,t,a,K);
                const mass = total_mass(R,L,b,t,a);
                const sf_loc = local_sf(R,L,b,t,a,N,K);
                const sf_glob = (N && N>0) ? (ncr/Math.max(N,1e-12)) : Infinity;
                const sf_min_val = Math.min(sf_loc, sf_glob);
                const ncr_over_m = ncr/Math.max(mass,1e-12);
                let objVal = (obj==='min_mass') ? mass : -ncr_over_m;
                if (sf_min_val < sfReq) objVal += 1e9;
                if (objVal < bestObj){ bestObj=objVal; best={xi:Xs[i], yi:Ys[j], zi:Zs[k], ncr, mass, sf_min:sf_min_val, ncr_over_m, b,t,a,L,nth, material:mname, lattice:lat}; }
              }
            }
          }
        state.lattice=latSaved; updateGeom(state);
      }
      mat = matSaved; state.material = stateMatSaved;
    }
    if(best){
      // Update UI to best selections
      state.material = best.material; mat = MAT_DB[state.material]; document.getElementById('material').value = state.material;
      state.lattice = best.lattice; document.getElementById('lattice').value = state.lattice; updateGeom(state);
      updateGeom(state);
      // Marker position
      const ax=state.xAxis, ay=state.yAxis;
      const R=mm_to_m(state.R_mm), K=state.K, KDF=state.KDF, N=state.N_req, M=state.M_req, T=state.T_req;
      const metricVal = METRICS[state.metric].fn(R, best.L, best.b, best.t, best.a, K, KDF, N, M, T);
      state.bestMarker = { x:(ax==='b'||ax==='t'||ax==='a')? best.xi*1000.0:best.xi, y:(ay==='b'||ay==='t'||ay==='a')? best.yi*1000.0:best.yi, z:metricVal };
      const li = (s)=>`<li>${s}</li>`;
      const AXL=AX_LABEL; const fmt=(name,val)=> (name==='b'||name==='t'||name==='a')? `${(val*1000).toFixed(2)} mm` : (name==='n_theta'? `${val}` : `${val.toFixed(3)}`);
      const bx=fmt(state.xAxis, best.xi), by=fmt(state.yAxis, best.yi), bz=fmt(state.zAxis, best.zi);
      const line1 = `Position: x=${AXL[state.xAxis]}=${bx}; y=${AXL[state.yAxis]}=${by}; Slice(${AXL[state.zAxis]})=${bz}`;
      const bmm=(best.b*1000).toFixed(2), tmm=(best.t*1000).toFixed(2), amm=(best.a*1000).toFixed(2);
      const line2 = `Geometrie: b=${bmm} mm; t=${tmm} mm; a=${amm} mm; n_theta=${best.nth}; L=${best.L.toFixed(3)} m; lattice=${best.lattice}; material=${best.material}`;
      const line3 = `Metriken: mass=${best.mass.toFixed(3)} kg; Ncr=${best.ncr.toExponential(2)} N; Ncr/m=${best.ncr_over_m.toFixed(3)} m/s²; SF_min=${best.sf_min.toFixed(2)}`;
      document.getElementById('optOut').innerHTML = `<ul style=\"margin:6px 0 0 18px;\">${li(line1)}${li(line2)}${li(line3)}</ul>`;
      rerender();
    } else {
      document.getElementById('optOut').textContent = 'Keine zulässige Kombination gefunden.';
    }
  }
}

window.addEventListener('load', mount);
</script>
</body>
</html>
"""


def main():
    ap = argparse.ArgumentParser(description="Version 2 Isogrid/Honeycomb Dashboard generator")
    ap.add_argument("--out", type=Path, default=Path("out/isogrid_dashboard_v2.html"))
    ap.add_argument("--R_mm", type=float, default=550.0)
    ap.add_argument("--K", type=float, default=1.0)
    ap.add_argument("--N_req", type=float, default=150000.0)
    ap.add_argument("--M_req", type=float, default=500.0)
    ap.add_argument("--T_req", type=float, default=800.0)
    ap.add_argument("--b_mm", nargs=2, type=float, default=[3.0, 15.0])
    ap.add_argument("--t_mm", nargs=2, type=float, default=[0.5, 4.0])
    ap.add_argument("--n_theta", nargs=2, type=int, default=[40, 140])
    ap.add_argument("--a_mm", nargs=2, type=float, default=[20.0, 70.0])
    ap.add_argument("--L", nargs=2, type=float, default=[0.3, 1.2])
    ap.add_argument("--nx", type=int, default=15)
    ap.add_argument("--ny", type=int, default=12)
    ap.add_argument("--nz", type=int, default=12)
    # lattice selection via dropdown; orientation handled in JS; no theta arg
    args = ap.parse_args()

    state = {
        "xAxis":"b", "yAxis":"t", "zAxis":"n_theta",
        "plotType":"surface",
        "material":"AISI 304",
        "lattice":"tri_0",
        "metric":"ncr_global",
        "KDF":0.85,
        "nx":args.nx, "ny":args.ny, "nz":args.nz,
        "b_min_mm":args.b_mm[0], "b_max_mm":args.b_mm[1],
        "t_min_mm":args.t_mm[0], "t_max_mm":args.t_mm[1],
        "n_min":args.n_theta[0], "n_max":args.n_theta[1],
        "a_min_mm":args.a_mm[0], "a_max_mm":args.a_mm[1],
        "L_min":args.L[0], "L_max":args.L[1],
        "b_fix_mm": (args.b_mm[0]+args.b_mm[1])/2,
        "t_fix_mm": (args.t_mm[0]+args.t_mm[1])/2,
        "n_fix": (args.n_theta[0]+args.n_theta[1])/2,
        "a_fix_mm": (args.a_mm[0]+args.a_mm[1])/2,
        "L_fix": (args.L[0]+args.L[1])/2,
        "R_mm": args.R_mm, "K": args.K,
        "N_req": args.N_req, "M_req": args.M_req, "T_req": args.T_req,
    }

    html = HTML_TEMPLATE % {
        "nx": args.nx, "ny": args.ny, "nz": args.nz,
        "b_min_mm": args.b_mm[0], "b_max_mm": args.b_mm[1],
        "t_min_mm": args.t_mm[0], "t_max_mm": args.t_mm[1],
        "n_min": args.n_theta[0], "n_max": args.n_theta[1],
        "a_min_mm": args.a_mm[0], "a_max_mm": args.a_mm[1],
        "L_min": args.L[0], "L_max": args.L[1],
        "R_mm": args.R_mm, "K": args.K, "KDF": 0.85,
        "N_req": args.N_req, "M_req": args.M_req, "T_req": args.T_req,
        "state_json": json.dumps(state),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html, encoding='utf-8')
    print(f"Wrote dashboard v2: {args.out}")


if __name__ == '__main__':
    main()
