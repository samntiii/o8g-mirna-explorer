#!/usr/bin/env bash
# Back-compat wrapper — prefer ./start.sh
exec "$(dirname "$0")/start.sh" "$@"
