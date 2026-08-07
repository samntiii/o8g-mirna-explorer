#!/usr/bin/env bash
# Interactive one-time Cloudflare Tunnel setup for samnti.com
set -euo pipefail
TUNNEL_NAME="${TUNNEL_NAME:-o8g-explorer}"
HOSTNAME="${HOSTNAME_PUBLIC:-oxomir.samnti.com}"

if ! command -v cloudflared >/dev/null 2>&1; then
  echo "Installing cloudflared via Homebrew…"
  brew install cloudflared
fi

echo ">>> Browser will open. Log into the Cloudflare account that owns $HOSTNAME."
cloudflared tunnel login

if cloudflared tunnel list 2>/dev/null | grep -q "$TUNNEL_NAME"; then
  echo "Tunnel '$TUNNEL_NAME' already exists."
else
  cloudflared tunnel create "$TUNNEL_NAME"
fi

echo ">>> Routing DNS (idempotent if already routed)…"
cloudflared tunnel route dns "$TUNNEL_NAME" "$HOSTNAME" || true
cloudflared tunnel route dns "$TUNNEL_NAME" "www.$HOSTNAME" || true

UUID=$(cloudflared tunnel list | awk -v n="$TUNNEL_NAME" '$0 ~ n {print $1; exit}')
CREDS="$HOME/.cloudflared/${UUID}.json"
CFG="$HOME/.cloudflared/config.yml"

mkdir -p "$HOME/.cloudflared"
cat > "$CFG" <<EOF
tunnel: ${UUID}
credentials-file: ${CREDS}

ingress:
  - hostname: ${HOSTNAME}
    service: http://127.0.0.1:8501
  - hostname: www.${HOSTNAME}
    service: http://127.0.0.1:8501
  - service: http_status:404
EOF

echo
echo "Wrote $CFG"
echo "Tunnel UUID: $UUID"
echo "Next:  cd ~/o8g_mirna_explorer && ./start_public.sh"
echo "Site:  https://${HOSTNAME}"
