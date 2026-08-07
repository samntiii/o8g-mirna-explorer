#!/usr/bin/env bash
# Public launch: Streamlit + Cloudflare Tunnel (oxomir.samnti.com).
# Children are detached — they keep running after this script exits.
# Stop with: stop_oxomir
set -euo pipefail
cd "$(dirname "$0")"

TUNNEL_NAME="${TUNNEL_NAME:-o8g-explorer}"
PORT="${PORT:-8501}"
CF_CONFIG="${CLOUDFLARED_CONFIG:-$(pwd)/deploy/cloudflared.yml}"
PUBLIC_HOST="${PUBLIC_HOST:-oxomir.samnti.com}"

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

mkdir -p pids

# Replace any prior instance
pkill -f 'cloudflared tunnel.*o8g-explorer' 2>/dev/null || true
pkill -f 'streamlit run app.py' 2>/dev/null || true
pkill -f 'caffeinate -dims -w' 2>/dev/null || true
if lsof -tiTCP:"$PORT" -sTCP:LISTEN >/dev/null 2>&1; then
  lsof -tiTCP:"$PORT" -sTCP:LISTEN | xargs kill -9 2>/dev/null || true
fi
sleep 0.8

echo "Starting Streamlit on :$PORT"
nohup "$ST" run app.py \
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
nohup cloudflared tunnel --config "$CF_CONFIG" run "$TUNNEL_NAME" \
  > pids/cloudflared.log 2>&1 &
echo $! > pids/cloudflared.pid
sleep 2
if ! pgrep -f "cloudflared tunnel.*$TUNNEL_NAME" >/dev/null 2>&1; then
  echo "ERROR: cloudflared failed to stay up" >&2
  tail -30 pids/cloudflared.log >&2 || true
  exit 1
fi

# Keep Mac awake while Streamlit lives (detached)
nohup caffeinate -dims -w "$(cat pids/streamlit.pid)" > pids/caffeinate.log 2>&1 &
echo $! > pids/caffeinate.pid

echo "Public site: https://$PUBLIC_HOST"
echo "Local:       http://127.0.0.1:$PORT"
echo "Stop with:   stop_oxomir"
