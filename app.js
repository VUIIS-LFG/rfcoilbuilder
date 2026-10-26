// app.js - builds the forms, sends settings to the Python worker (worker.js) and shows the results.

// ---------------------------------------------------------------- form definitions
const PROFILE_OPTIONS = [
  ['uniform', 'Uniform (flat)'], ['sqrt', 'Square root (linear gradient)'], ['linear', 'Linear gradient'],
  ['quadratic', 'Quadratic'], ['cosine', 'Cosine bell'], ['gaussian', 'Gaussian'], ['linear_symmetric', 'Linear symmetric |x|'],
];
const SPACING_OPTIONS = [['uniform', 'Uniform'], ['sqrt', 'Square root (denser at the ends)'], ['cosine', 'Cosine (denser at the ends)']];

const FORMS = {
  solenoid: [
    { group: 'Coil shape', fields: [
      { id: 'coil_length_mm', label: 'Coil length', unit: 'mm', v: 160, min: 10 },
      { id: 'coil_radius_mm', label: 'Inner radius', unit: 'mm', v: 50, min: 5, help: 'The inner diameter is twice this.' },
      { id: 'former_thickness_mm', label: 'Former wall thickness', unit: 'mm', v: 3, min: 0.5 },
    ] },
    { group: 'Wire', fields: [
      { id: 'wire_bare_mm', label: 'Bare wire diameter', unit: 'mm', v: 0.812, min: 0.05, help: 'AWG 20 is 0.812 mm bare and 0.879 mm insulated. Check your wire\'s data sheet.' },
      { id: 'wire_ins_mm', label: 'Insulated wire diameter', unit: 'mm', v: 0.879, min: 0.05 },
    ] },
    { group: 'Grooves and turns', fields: [
      { id: 'M', label: 'Number of grooves', v: 24, min: 2, max: 60, step: 1 },
      { id: 'spacing', label: 'Groove spacing', type: 'select', options: SPACING_OPTIONS, v: 'uniform' },
      { id: 'n_max', label: 'Max forward turns per groove', v: 2, min: 1, max: 6, step: 1 },
      { id: 'n_min', label: 'Max bucking turns per groove', v: -2, min: -6, max: 0, step: 1, help: 'Enter as a negative number.' },
      { id: 'bucking', label: 'Allow bucking (reverse) turns', type: 'check', v: true },
    ] },
    { group: 'Target field', fields: [
      { id: 'profile', label: 'Field profile along the axis', type: 'select', options: PROFILE_OPTIONS, v: 'uniform' },
      { id: 'sigma', label: 'Gaussian width (only for Gaussian)', v: 0.25, min: 0.05, max: 1, step: 0.05 },
    ] },
    { group: 'Operation', fields: [
      { id: 'freq_mhz', label: 'Operating frequency', unit: 'MHz', v: 2.075, min: 0.01 },
    ] },
    { group: 'Former (used for the build instructions)', fields: [
      { id: 'end_margin_mm', label: 'Extra former beyond the end grooves', unit: 'mm', v: 10, min: 0 },
      { id: 'groove_depth_mm', label: 'Groove depth', unit: 'mm', v: 1.0, min: 0.1, help: 'Must be less than the wall thickness.' },
      { id: 'groove_clearance_mm', label: 'Groove clearance', unit: 'mm', v: 0.2, min: 0, help: 'Groove width = insulated wire diameter + this.' },
    ] },
    { group: 'Calculation', open: false, fields: [
      { id: 'N', label: 'Field points along the axis', v: 60, min: 10, max: 400, step: 1 },
      { id: 'map_res_mm', label: 'Field map resolution', unit: 'mm', v: 2, min: 1, help: 'Smaller is sharper but slower.' },
      { id: 'solver_s', label: 'Optimiser time limit', unit: 's', v: 20, min: 2, max: 300, step: 1 },
    ] },
  ],
  conformal: [
    { group: 'Phantom', fields: [
      { id: 'stl', label: 'Phantom STL file', type: 'file', help: 'Leave empty to try the built-in test head. Your file stays in your browser and is not uploaded anywhere.' },
      { id: 'units_to_m', label: 'File units', type: 'select', options: [['0.001', 'Millimetres'], ['0.01', 'Centimetres'], ['1', 'Metres'], ['0.0254', 'Inches']], v: '0.001' },
      { id: 'axis', label: 'Axis that runs along the head', type: 'select', options: [['auto', 'Automatic (longest dimension)'], ['x', 'X'], ['y', 'Y'], ['z', 'Z']], v: 'auto', help: 'This becomes the coil axis.' },
    ] },
    { group: 'Former (shell)', fields: [
      { id: 'standoff_mm', label: 'Gap between head and former', unit: 'mm', v: 10, min: 0 },
      { id: 'former_thickness_mm', label: 'Former wall thickness', unit: 'mm', v: 3, min: 0.5 },
      { id: 'slide_from', label: 'The coil slides onto the head from', type: 'select', options: [['top', 'The top'], ['bottom', 'The neck side'], ['none', 'Neither (split former)']], v: 'top', help: 'Widens the shell where needed so the head can pass through.' },
      { id: 'groove_depth_mm', label: 'Groove depth', unit: 'mm', v: 1.0, min: 0.1 },
      { id: 'groove_clearance_mm', label: 'Groove clearance', unit: 'mm', v: 0.2, min: 0 },
      { id: 'end_margin_mm', label: 'Extra former beyond the end grooves', unit: 'mm', v: 10, min: 0 },
    ] },
    { group: 'Wire', fields: [
      { id: 'wire_bare_mm', label: 'Bare wire diameter', unit: 'mm', v: 0.812, min: 0.05 },
      { id: 'wire_ins_mm', label: 'Insulated wire diameter', unit: 'mm', v: 0.879, min: 0.05 },
    ] },
    { group: 'Grooves and turns', fields: [
      { id: 'M', label: 'Number of grooves', v: 24, min: 2, max: 60, step: 1 },
      { id: 'spacing', label: 'Groove spacing', type: 'select', options: SPACING_OPTIONS, v: 'uniform' },
      { id: 'n_max', label: 'Max forward turns per groove', v: 2, min: 1, max: 6, step: 1 },
      { id: 'n_min', label: 'Max bucking turns per groove', v: -2, min: -6, max: 0, step: 1 },
      { id: 'bucking', label: 'Allow bucking (reverse) turns', type: 'check', v: false },
      { id: 'target_fraction', label: 'Target field strength', v: 1.0, min: 0.1, max: 3, step: 0.1, help: '1.0 = the field a one-turn-per-groove cylinder would give.' },
      { id: 'freq_mhz', label: 'Operating frequency', unit: 'MHz', v: 2.075, min: 0.01 },
    ] },
    { group: 'Imaging region', fields: [
      { id: 'fov_mode', label: 'Region the field should be uniform over', type: 'select', options: [['anatomy', 'The head\'s own shape'], ['box', 'A plain box']], v: 'anatomy' },
      { id: 'fov_x_mm', label: 'Region size x', unit: 'mm', v: '', placeholder: 'auto', min: 1, help: 'Leave x, y and z empty for automatic sizing (60%, 60% and 40% of the head).' },
      { id: 'fov_y_mm', label: 'Region size y', unit: 'mm', v: '', placeholder: 'auto', min: 1 },
      { id: 'fov_z_mm', label: 'Region size z (along the coil)', unit: 'mm', v: '', placeholder: 'auto', min: 1 },
    ] },
    { group: 'Advanced', open: false, fields: [
      { id: 'fov_cx_mm', label: 'Region centre x', unit: 'mm', v: 0 },
      { id: 'fov_cy_mm', label: 'Region centre y', unit: 'mm', v: 0 },
      { id: 'fov_cz_mm', label: 'Region centre z', unit: 'mm', v: 0 },
      { id: 'coil_margin_mm', label: 'Coil margin beyond the region (each end)', unit: 'mm', v: '', placeholder: 'auto', min: 0 },
      { id: 'roi_shrink', label: 'Fraction of the head radius to use', v: 0.85, min: 0.3, max: 1, step: 0.05 },
      { id: 'roi_spacing_mm', label: 'Sample point spacing', unit: 'mm', v: 15, min: 4, help: 'Smaller means more sample points and a slower design.' },
      { id: 'solver_s', label: 'Optimiser time limit (per coil)', unit: 's', v: 20, min: 2, max: 300, step: 1 },
    ] },
  ],
};

// ---------------------------------------------------------------- helpers
const $ = (sel, root = document) => root.querySelector(sel);
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

function renderForm(which) {
  const form = $(`#form-${which}`);
  form.innerHTML = FORMS[which].map((g) => {
    const open = g.open === false ? '' : ' open';
    const rows = g.fields.map((f) => {
      const id = `${which}-${f.id}`;
      const help = f.help ? `<span class="help" id="${id}-help">${esc(f.help)}</span>` : '';
      const desc = f.help ? ` aria-describedby="${id}-help"` : '';
      if (f.type === 'select') {
        const opts = f.options.map(([val, text]) => `<option value="${val}"${val === f.v ? ' selected' : ''}>${esc(text)}</option>`).join('');
        return `<div class="field"><label for="${id}">${esc(f.label)}</label><select id="${id}"${desc}>${opts}</select>${help}</div>`;
      }
      if (f.type === 'check') {
        return `<div class="field check"><input type="checkbox" id="${id}"${f.v ? ' checked' : ''}${desc}><label for="${id}">${esc(f.label)}</label>${help}</div>`;
      }
      if (f.type === 'file') {
        return `<div class="field"><label for="${id}">${esc(f.label)}</label><input type="file" id="${id}" accept=".stl"${desc}>${help}</div>`;
      }
      const attrs = [`min="${f.min ?? ''}"`, f.max !== undefined ? `max="${f.max}"` : '', `step="${f.step ?? 'any'}"`,
        f.placeholder ? `placeholder="${f.placeholder}"` : ''].join(' ');
      const unit = f.unit ? `<span class="unit">${esc(f.unit)}</span>` : '';
      return `<div class="field"><label for="${id}">${esc(f.label)}</label><div class="num"><input type="number" id="${id}" value="${f.v}" ${attrs}${desc}>${unit}</div>${help}</div>`;
    }).join('');
    return `<details class="group"${open}><summary>${esc(g.group)}</summary>${rows}</details>`;
  }).join('');
}

function collect(which) {
  const params = {};
  for (const g of FORMS[which]) {
    for (const f of g.fields) {
      if (f.type === 'file') continue;
      const el = $(`#${which}-${f.id}`);
      if (f.type === 'check') params[f.id] = el.checked;
      else if (f.type === 'select') params[f.id] = el.value;
      else params[f.id] = el.value === '' ? '' : Number(el.value);
    }
  }
  if (which === 'conformal') params.units_to_m = Number(params.units_to_m);
  return params;
}

// ---------------------------------------------------------------- worker
const worker = new Worker('worker.js');
let engineReady = false;
let jobId = 0;
const lastResult = {};

worker.onmessage = (e) => {
  const m = e.data;
  if (m.type === 'status') setEngine(m.text, false);
  else if (m.type === 'ready') {
    engineReady = true;
    setEngine('Python engine ready', true);
    document.querySelectorAll('.run').forEach((b) => (b.disabled = false));
  } else if (m.type === 'initerror') {
    setEngine('The Python engine failed to load: ' + m.message, false, true);
  } else if (m.type === 'result' && m.id === jobId) {
    finish(m.which);
    if (m.data.error) showError(m.which, m.data.error);
    else { lastResult[m.which] = m.data; showResults(m.which, m.data); }
  } else if (m.type === 'error' && m.id === jobId) {
    finish(m.which);
    showError(m.which, m.message);
  }
};
worker.onerror = (e) => setEngine('The Python engine hit an error: ' + (e.message || 'unknown'), false, true);

function setEngine(text, ok, bad) {
  const el = $('#engine');
  el.textContent = text;
  el.className = 'engine' + (ok ? ' ok' : '') + (bad ? ' bad' : '');
}

function finish(which) {
  $(`#run-${which}`).disabled = false;
  $(`#status-${which}`).textContent = '';
  $(`#panel-${which}`).setAttribute('aria-busy', 'false');
}

async function run(which) {
  if (!engineReady) return;
  const form = $(`#form-${which}`);
  const bad = [...form.querySelectorAll('input[type=number]')].find((i) => !i.checkValidity());
  if (bad) { bad.reportValidity(); return; }
  const params = collect(which);
  let file, transfer = [];
  if (which === 'conformal') {
    const f = $('#conformal-stl').files[0];
    if (f) { file = await f.arrayBuffer(); transfer = [file]; }
  }
  $(`#run-${which}`).disabled = true;
  $(`#panel-${which}`).setAttribute('aria-busy', 'true');
  $(`#status-${which}`).textContent = which === 'conformal'
    ? 'Designing the coil. This usually takes under a minute, longer for large STL files.'
    : 'Designing the coil. This usually takes about half a minute.';
  $(`#results-${which}`).innerHTML = '';
  jobId += 1;
  worker.postMessage({ id: jobId, which, params, file }, transfer);
}

// ---------------------------------------------------------------- results
function table(head, rows, rowClass) {
  const th = head.map((h) => `<th scope="col">${esc(h)}</th>`).join('');
  const tr = rows.map((r, i) => `<tr${rowClass ? ` class="${rowClass(i)}"` : ''}>${r.map((c) => `<td>${esc(c)}</td>`).join('')}</tr>`).join('');
  return `<div class="scroll"><table><thead><tr>${th}</tr></thead><tbody>${tr}</tbody></table></div>`;
}

function showError(which, message) {
  const line = String(message).trim().split('\n').pop();
  $(`#results-${which}`).innerHTML = `<div class="callout bad" role="alert"><strong>The design did not run.</strong> ${esc(line)}</div>`;
}

function showResults(which, d) {
  let html = '';
  if (d.notes && d.notes.length) {
    html += `<div class="callout" role="status"><strong>Worth knowing</strong><ul>${d.notes.map((n) => `<li>${esc(n)}</li>`).join('')}</ul></div>`;
  }
  html += `<section><h3>Summary</h3><dl class="summary">${d.summary.map(([k, v]) => `<div><dt>${esc(k)}</dt><dd>${esc(v)}</dd></div>`).join('')}</dl></section>`;
  if (d.comparison) {
    html += `<section><h3>Conformal coil vs cylinder</h3>${table(['', 'Conformal', 'Cylinder'], d.comparison)}</section>`;
  }
  // winding pattern
  if (which === 'solenoid') {
    html += `<section><h3>Winding pattern</h3>${table(['Groove', 'From former end (mm)', 'From coil centre (mm)', 'Turns', 'Direction'],
      d.windings.map((w) => [w.groove, w.from_end_mm.toFixed(2), w.from_center_mm.toFixed(2), (w.turns > 0 ? '+' : '') + w.turns, w.direction]),
      (i) => (d.windings[i].turns > 0 ? 'fwd' : 'bck'))}`;
  } else {
    html += `<section><h3>Winding pattern</h3>${table(['Groove', 'From former end (mm)', 'Conformal turns', 'Cylinder turns'],
      d.windings.map((w) => [w.groove, w.from_end_mm.toFixed(2), sign(w.conformal), sign(w.cylinder)]))}`;
  }
  html += `<p><button type="button" class="secondary" data-csv="${which}">Download winding pattern (CSV)</button></p></section>`;
  html += `<section><h3>How to build the former</h3><ol class="steps">${d.instructions.map((s) => `<li>${esc(s)}</li>`).join('')}</ol></section>`;
  html += `<section><h3>Graphs</h3>${d.figures.map((f, i) => `<figure><img alt="${esc(f.title + '. ' + f.caption)}" src="data:image/png;base64,${f.png}">
    <figcaption><strong>${esc(f.title)}.</strong> ${esc(f.caption)} <a download="${which}-${i + 1}-${f.title.toLowerCase().replace(/[^a-z0-9]+/g, '-')}.png" href="data:image/png;base64,${f.png}">Download PNG</a></figcaption></figure>`).join('')}</section>`;
  $(`#results-${which}`).innerHTML = html;
  $(`#results-${which}`).scrollIntoView({ behavior: 'smooth', block: 'start' });
}

const sign = (n) => (n > 0 ? '+' + n : String(n));

function downloadCsv(which) {
  const d = lastResult[which];
  if (!d) return;
  const rows = which === 'solenoid'
    ? [['groove', 'from_former_end_mm', 'from_coil_center_mm', 'turns', 'direction'], ...d.windings.map((w) => [w.groove, w.from_end_mm, w.from_center_mm, w.turns, w.direction])]
    : [['groove', 'from_former_end_mm', 'conformal_turns', 'cylinder_turns'], ...d.windings.map((w) => [w.groove, w.from_end_mm, w.conformal, w.cylinder])];
  const blob = new Blob([rows.map((r) => r.join(',')).join('\n')], { type: 'text/csv' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `${which}-winding-pattern.csv`;
  a.click();
  URL.revokeObjectURL(a.href);
}

// ---------------------------------------------------------------- tabs, buttons, header motif
function selectTab(name) {
  document.querySelectorAll('[role=tab]').forEach((t) => {
    const on = t.dataset.tab === name;
    t.setAttribute('aria-selected', on);
    t.tabIndex = on ? 0 : -1;
    $(`#panel-${t.dataset.tab}`).hidden = !on;
  });
  history.replaceState(null, '', '#' + name);
}

function init() {
  renderForm('solenoid');
  renderForm('conformal');
  const tabs = [...document.querySelectorAll('[role=tab]')];
  tabs.forEach((t, i) => {
    t.addEventListener('click', () => selectTab(t.dataset.tab));
    t.addEventListener('keydown', (e) => {
      const j = e.key === 'ArrowRight' ? i + 1 : e.key === 'ArrowLeft' ? i - 1 : null;
      if (j === null) return;
      const next = tabs[(j + tabs.length) % tabs.length];
      next.focus(); selectTab(next.dataset.tab);
    });
  });
  const start = location.hash.slice(1);
  selectTab(['solenoid', 'conformal', 'development'].includes(start) ? start : 'solenoid');

  document.querySelectorAll('.run').forEach((b) => { b.disabled = true; b.addEventListener('click', () => run(b.dataset.which)); });
  document.querySelectorAll('.reset').forEach((b) => b.addEventListener('click', () => {
    renderForm(b.dataset.which);
    $(`#results-${b.dataset.which}`).innerHTML = '';
  }));
  document.addEventListener('click', (e) => { const b = e.target.closest('[data-csv]'); if (b) downloadCsv(b.dataset.csv); });

  // header motif: wire cross-sections above and below a tube, in the style of the side view graph
  const turns = [1, 2, 2, -1, 2, 1, -1, 0, -1, 1, 2, -1, 2, 2, 1];
  let svg = '<line x1="0" y1="30" x2="300" y2="30" stroke="#7a8793" stroke-dasharray="6 4"/><rect x="6" y="19" width="288" height="3" fill="#E69F00"/><rect x="6" y="38" width="288" height="3" fill="#E69F00"/>';
  turns.forEach((n, i) => {
    const x = 14 + i * 19.5;
    const c = n > 0 ? '#0072B2' : '#D55E00';
    for (let k = 0; k < Math.abs(n); k++) {
      svg += `<circle cx="${x}" cy="${16 - k * 6}" r="3.2" fill="${c}"/><circle cx="${x}" cy="${44 + k * 6}" r="3.2" fill="${c}"/>`;
    }
  });
  $('#motif').innerHTML = svg;
}
init();
