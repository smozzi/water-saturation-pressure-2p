#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "This script must be sourced so the virtualenv stays active in your shell." >&2
    echo "Usage: source scripts/venv.sh" >&2
    exit 1
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
venv_dir="${repo_root}/.venv"

if [[ ! -d "${venv_dir}" ]]; then
    echo "Creating virtual environment in ${venv_dir}" >&2
    python -m venv "${venv_dir}"
fi

# shellcheck disable=SC1090
source "${venv_dir}/bin/activate"

echo "Activated virtual environment at ${venv_dir}" >&2
