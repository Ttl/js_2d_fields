// Resizable panes of the main column: the log panel under the tabs and the custom
// geometry editor under the preview. Sizes and the log state persist in localStorage.

const store = {
    get(key) { try { return localStorage.getItem(key); } catch { return null; } },
    set(key, value) { try { localStorage.setItem(key, value); } catch { /* unavailable */ } },
};

const LOG_COLLAPSED_KEY = 'log_collapsed';
const LOG_HEIGHT_KEY = 'log_height';
const EDITOR_SHARE_KEY = 'custom_editor_share';

const LOG_MIN = 40;
const PREVIEW_MIN = 120;
const EDITOR_MIN = 140;

let resizeQueued = false;
function notifyResize() {
    if (resizeQueued) return;
    resizeQueued = true;
    requestAnimationFrame(() => { resizeQueued = false; window.dispatchEvent(new Event('resize')); });
}

// Vertical drag on a handle. onStart runs at the grab, onMove gets the upward travel
// since then in px, onEnd runs once at release.
function dragHandle(handle, onStart, onMove, onEnd) {
    handle.addEventListener('pointerdown', (e) => {
        if (e.button !== undefined && e.button !== 0) return;
        e.preventDefault();
        handle.setPointerCapture(e.pointerId);
        handle.classList.add('dragging');
        const y0 = e.clientY;
        onStart();
        const move = (ev) => { onMove(y0 - ev.clientY); notifyResize(); };
        const up = () => {
            handle.classList.remove('dragging');
            handle.removeEventListener('pointermove', move);
            handle.removeEventListener('pointerup', up);
            handle.removeEventListener('pointercancel', up);
            if (onEnd) onEnd();
            notifyResize();
        };
        handle.addEventListener('pointermove', move);
        handle.addEventListener('pointerup', up);
        handle.addEventListener('pointercancel', up);
    });
}

// --- Log panel ----------------------------------------------------------------------

const isCustomType = () => document.getElementById('tl_type')?.value === 'custom';

// The user's own choice wins. Without one the log is collapsed while a custom geometry
// is being edited and opens when a solve starts.
let autoOpened = false;

function logCollapsed() {
    const saved = store.get(LOG_COLLAPSED_KEY);
    if (saved !== null) return saved === '1';
    return isCustomType() && !autoOpened;
}

export function syncLogPanel() {
    const panel = document.getElementById('log-panel');
    if (!panel) return;
    const collapsed = logCollapsed();
    if (panel.classList.contains('collapsed') === collapsed) return;
    panel.classList.toggle('collapsed', collapsed);
    const btn = document.getElementById('btn-log-toggle');
    btn.setAttribute('aria-expanded', String(!collapsed));
    btn.title = collapsed ? 'Show the log' : 'Hide the log';
    if (!collapsed) {
        const out = document.getElementById('console_out');
        out.scrollTop = out.scrollHeight;
    }
    notifyResize();
}

// Last log line, shown in the bar while the log is collapsed.
export function setLogStatus(line) {
    const status = document.getElementById('log-last');
    if (!status) return;
    status.textContent = line;
    status.classList.toggle('is-error', /^(ERROR|Error|Sweep error)/.test(line));
    status.classList.toggle('is-warning', line.startsWith('⚠'));
}

export function logSolveStarted() {
    autoOpened = true;
    syncLogPanel();
}

function initLogPanel() {
    const panel = document.getElementById('log-panel');
    if (!panel) return;
    const height = parseFloat(store.get(LOG_HEIGHT_KEY));
    if (Number.isFinite(height)) panel.style.setProperty('--log-height', `${height}px`);

    document.getElementById('log-bar').addEventListener('click', () => {
        store.set(LOG_COLLAPSED_KEY, logCollapsed() ? '0' : '1');
        syncLogPanel();
    });
    document.getElementById('tl_type').addEventListener('change', () => { autoOpened = false; syncLogPanel(); });

    const out = document.getElementById('console_out');
    let start = 0, current = 0;
    dragHandle(document.getElementById('log-splitter'), () => { start = out.offsetHeight; }, (up) => {
        if (panel.classList.contains('collapsed')) return;
        const max = panel.parentElement.clientHeight * 0.7;
        current = Math.min(max, Math.max(LOG_MIN, start + up));
        panel.style.setProperty('--log-height', `${current}px`);
    }, () => { if (current) store.set(LOG_HEIGHT_KEY, String(Math.round(current))); });

    syncLogPanel();
}

// --- Custom geometry editor ---------------------------------------------------------

function initEditorSplitter() {
    const editor = document.getElementById('custom-editor');
    const handle = document.getElementById('custom-splitter');
    const tab = document.getElementById('tab-geometry');
    if (!editor || !handle || !tab) return;

    const setShare = (share) => { editor.style.setProperty('--editor-share', `${(share * 100).toFixed(1)}%`); };
    const saved = parseFloat(store.get(EDITOR_SHARE_KEY));
    if (saved > 0 && saved < 1) setShare(saved);

    let start = 0, share = 0;
    dragHandle(handle, () => { start = editor.offsetHeight; }, (up) => {
        if (editor.classList.contains('collapsed')) return;
        const total = tab.clientHeight;
        const px = Math.min(total - PREVIEW_MIN, Math.max(EDITOR_MIN, start + up));
        share = px / total;
        setShare(share);
    }, () => { if (share) store.set(EDITOR_SHARE_KEY, share.toFixed(3)); });
}

export function initLayoutPanels() {
    initLogPanel();
    initEditorSplitter();
}
