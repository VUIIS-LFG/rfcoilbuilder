// worker.js - runs Python in the browser (Pyodide) so the page stays responsive while a coil is designed.
importScripts('https://cdn.jsdelivr.net/pyodide/v0.26.4/full/pyodide.js');

const engines = {};

async function init() {
  postMessage({ type: 'status', text: 'Loading the Python engine (first visit can take a minute or two)...' });
  const pyodide = await loadPyodide();
  postMessage({ type: 'status', text: 'Loading numpy, scipy and matplotlib...' });
  await pyodide.loadPackage(['numpy', 'scipy', 'matplotlib']);
  for (const name of ['solenoid', 'conformal']) {
    const response = await fetch(name + '.py');
    if (!response.ok) throw new Error('Could not load ' + name + '.py - is it uploaded next to index.html?');
    const source = await response.text();
    const scope = pyodide.globals.get('dict')();
    pyodide.runPython(source, { globals: scope });
    engines[name] = scope.get('run');
  }
  postMessage({ type: 'ready' });
  return pyodide;
}

const ready = init().catch((err) => {
  postMessage({ type: 'initerror', message: String(err) });
  throw err;
});
ready.catch(() => {});

self.onmessage = async (event) => {
  const { id, which, params, file } = event.data;
  try {
    const pyodide = await ready;
    let path = '';
    if (file) {
      pyodide.FS.writeFile('/phantom.stl', new Uint8Array(file));
      path = '/phantom.stl';
    }
    const out = engines[which](JSON.stringify(params), path);
    postMessage({ type: 'result', id, which, data: JSON.parse(out) });
  } catch (err) {
    postMessage({ type: 'error', id, which, message: String(err) });
  }
};
