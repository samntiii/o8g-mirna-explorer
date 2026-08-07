# Deploying o8G-miRNA Explorer to samnti.com

The app is a Streamlit server with an ~90 MB local SQLite database. A public
URL needs an always-reachable process. Two routes:

1. **Cloudflare Tunnel from this Mac** (current setup) — free, custom domain,
   Mac must stay awake.
2. **Render / Railway / Fly** via the `Dockerfile` — always-on, ~$5–7/mo.

Shiny / shinyapps.io is **not** simpler for a tunnel: `cloudflared` forwards
any local HTTP server. Rewriting the Python engine in R would not help hosting.

---

## Option A — Cloudflare Tunnel → https://oxomir.samnti.com

### One-time (interactive)

1. Put **samnti.com** on Cloudflare DNS (nameservers at the registrar).
2. Install the agent:

```bash
brew install cloudflared
```

3. Log in (opens a browser; pick the Cloudflare account that owns samnti.com):

```bash
cloudflared tunnel login
cloudflared tunnel create o8g-explorer
cloudflared tunnel route dns o8g-explorer oxomir.samnti.com
cloudflared tunnel route dns o8g-explorer o8g.samnti.com   # alias
cloudflared tunnel list    # copy the tunnel UUID
```

   **Note:** the apex `samnti.com` already has a DNS record (Foundry tunnel).
   This explorer is published at **https://oxomir.samnti.com** so the root site
   is left alone. To put the explorer on the apex later, delete/replace that
   A/CNAME in the Cloudflare dashboard, then `cloudflared tunnel route dns
   o8g-explorer samnti.com`.

4. Ingress config lives at
   [`deploy/cloudflared.yml`](deploy/cloudflared.yml) (already filled in on this
   machine). Template:
   [`deploy/cloudflared.config.example.yml`](deploy/cloudflared.config.example.yml).

### Run whenever the site should be live

```bash
cd ~/o8g_mirna_explorer
./start_public.sh
```

This starts Streamlit on `:8501`, runs the named tunnel, and `caffeinate`s the
Mac so it does not sleep. Ctrl-C stops both. Public URL:
**https://oxomir.samnti.com** (alias **https://o8g.samnti.com**).

Optional reboot persistence:

```bash
cp deploy/com.samnti.o8g-explorer.plist ~/Library/LaunchAgents/
launchctl load ~/Library/LaunchAgents/com.samnti.o8g-explorer.plist
```

Unload with `launchctl unload ~/Library/LaunchAgents/com.samnti.o8g-explorer.plist`.

**Caveats.** Sleep, Wi-Fi drops, quitting the script, or a laptop lid-close
without `caffeinate` take oxomir.samnti.com offline. Streamlit websockets are proxied
by cloudflared by default (`/_stcore/stream`).

---

## Option B — Render (always-on, no Mac)

1. Push this folder to GitHub (`o8g_targets.db` is ~87 MB, under the 100 MB
   GitHub file limit; use Git LFS if it grows).
2. Render → New → Web Service → connect the repo. It picks up the `Dockerfile`
   and `render.yaml`.
3. Custom domain: add `samnti.com` in Render, then at Cloudflare create the
   **CNAME / ALIAS** records Render shows (DNS only — not a tunnel).

Railway and Fly.io work the same way from the Dockerfile.

---

## Local only

```bash
cd ~/o8g_mirna_explorer
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./start.sh          # http://localhost:8501
```

Rebuild the database (only if inputs change):

```bash
python scripts/fetch_mirbase.py
python scripts/fetch_utr3.py          # skip if utr3_human.parquet exists
python precompute_db.py
python scripts/validate_db.py
```
