#!/usr/bin/env bash
# Public launch: Streamlit + Cloudflare Tunnel (oxomir.samnti.com).
#
# Default: stay in the foreground (wait on children) so launchd KeepAlive works.
# Fire-and-forget: DETACH=1 ./start_public.sh
# Stop: stop_oxomir   or   launchctl bootout gui/$(id -u)/com.samnti.o8g-explorer
set -euo pipefail
cd "$(dirname "$0")"

TUNNEL_NAME="${TUNNEL_NAME:-o8g-explorer}"
PORT="${PORT:-8501}"
CF_CONFIG="${CLOUDFLARED_CONFIG:-$(pwd)/deploy/cloudflared.yml}"
PUBLIC_HOST="${PUBLIC_HOST:-oxomir.samnti.com}"
DETACH="${DETACH:-0}"

export PATH="/opt/homebrew/bin:/usr/local/bin:${PATH:-}"

if ! command -v cloudflared >/dev/null 2>&1; then
  echo "ERROR: cloudflared not on PATH. Install with: brew install cloudflared" >&2
  exit 1
fi
if [ ! -f o8g_targets.db ]; then
  echo "ERROR: o8g_targets.db missing." >&2
  exit 1
fi
if [ ! -f "$CF_CONFIG" ]; then
  echo "ERROR: $CF_CONFIG not found. See DEPLOY.md." >&2
  exit 1
fi

if [ -x .venv/bin/streamlit ]; then
  ST=.venv/bin/streamlit
else
  ST=streamlit
fi

# ViennaRNA CLIs (RNAduplex / RNAup) for Gene→miRNA energetics
if [ -d /opt/homebrew/anaconda3/envs/mirna_viewer/bin ]; then
  export PATH="/opt/homebrew/anaconda3/envs/mirna_viewer/bin:$PATH"
fi

mkdir -p pids

cleanup() {
  local st_pid cf_pid caf_pid
  st_pid=$(cat pids/streamlit.pid 2>/dev/null || true)
  cf_pid=$(cat pids/cloudflared.pid 2>/dev/null || true)
  caf_pid=$(cat pids/caffeinate.pid 2>/dev/null || true)
  [ -n "${caf_pid:-}" ] && kill "$caf_pid" 2>/dev/null || true
  [ -n "${cf_pid:-}" ] && kill "$cf_pid" 2>/dev/null || true
  [ -n "${st_pid:-}" ] && kill "$st_pid" 2>/dev/null || true
  pkill -f "cloudflared tunnel.*${TUNNEL_NAME}" 2>/dev/null || true
  pkill -f 'streamlit run app.py' 2>/dev/null || true
}
trap cleanup EXIT INT TERM

# Replace any prior instance
pkill -f "cloudflared tunnel.*${TUNNEL_NAME}" 2>/dev/null || true
pkill -f 'streamlit run app.py' 2>/dev/null || true
pkill -f 'caffeinate -dims -w' 2>/dev/null || true
if lsof -tiTCP:"$PORT" -sTCP:LISTEN >/dev/null 2>&1; then
  lsof -tiTCP:"$PORT" -sTCP:LISTEN | xargs kill -9 2>/dev/null || true
fi
sleep 0.8

echo "Starting Streamlit on :$PORT"
"$ST" run app.py \
  --server.port="$PORT" \
  --server.address=0.0.0.0 \
  --server.headless=true \
  --server.enableCORS=false \
  --server.enableXsrfProtection=false \
  > pids/streamlit.log 2>&1 &
echo $! > pids/streamlit.pid

ready=0
for _ in $(seq 1 60); do
  if curl -sf "http://127.0.0.1:$PORT/_stcore/health" >/dev/null 2>&1; then
    ready=1
    break
  fi
  sleep 0.5
done
if [ "$ready" -ne 1 ]; then
  echo "ERROR: Streamlit did not become healthy on :$PORT" >&2
  tail -20 pids/streamlit.log >&2 || true
  exit 1
fi

echo "Starting cloudflared tunnel '$TUNNEL_NAME' → $PUBLIC_HOST"
cloudflared tunnel --config "$CF_CONFIG" run "$TUNNEL_NAME" \
  > pids/cloudflared.log 2>&1 &
echo $! > pids/cloudflared.pid
sleep 2
if ! pgrep -f "cloudflared tunnel.*$TUNNEL_NAME" >/dev/null 2>&1; then
  echo "ERROR: cloudflared failed to stay up" >&2
  tail -30 pids/cloudflared.log >&2 || true
  exit 1
fi

# Keep Mac awake while Streamlit lives
caffeinate -dims -w "$(cat pids/streamlit.pid)" > pids/caffeinate.log 2>&1 &
echo $! > pids/caffeinate.pid

echo "Public site: https://$PUBLIC_HOST"
echo "Local:       http://127.0.0.1:$PORT"
echo "Stop with:   stop_oxomir  (or unload LaunchAgent)"

if [ "$DETACH" = "1" ]; then
  trap - EXIT INT TERM
  echo "(detached — children keep running after this script exits)"
  exit 0
fi

# Stay alive so launchd KeepAlive does not restart-loop / kill children.
# Exit (and relaunch via KeepAlive) if either service dies.
st_pid=$(cat pids/streamlit.pid)
cf_pid=$(cat pids/cloudflared.pid)
while kill -0 "$st_pid" 2>/dev/null && kill -0 "$cf_pid" 2>/dev/null; do
  sleep 5
done
echo "A child process exited — shutting down so KeepAlive can restart." >&2
exit 1
