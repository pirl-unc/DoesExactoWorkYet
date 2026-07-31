#!/usr/bin/env bash
# Install the Exacto build under test and record exactly what it was.
#
# EXACTO_VERSION=latest-release (default) installs the newest published release
# tarball; anything else is treated as a git ref and built from source, which is
# how you point this repo at a branch or an unreleased commit.
set -euo pipefail

VERSION="${EXACTO_VERSION:-latest-release}"
REPO="${EXACTO_REPO:-https://github.com/pirl-unc/exacto}"
WORK="${DEWY_WORK_DIR:-work}/exacto-src"
RESULTS="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/results"

mkdir -p "$WORK" "$RESULTS"

api_repo="${REPO#https://github.com/}"

# api.github.com allows 60 unauthenticated requests per hour per IP, and CI
# runners share egress. At eight matrix legs that was invisible; at twenty-five
# starting together it returns 403 and every leg dies in its install step
# before doing any work. A token raises the limit to 5,000/hour and costs
# nothing -- GITHUB_TOKEN is already present in Actions.
gh_api() {
  local url="$1" attempt
  for attempt in 1 2 3; do
    if [[ -n "${GITHUB_TOKEN:-}" ]]; then
      curl -fsSL -H "Authorization: Bearer ${GITHUB_TOKEN}" \
           -H "X-GitHub-Api-Version: 2022-11-28" "$url" && return 0
    else
      curl -fsSL "$url" && return 0
    fi
    echo "  ${url} failed (attempt ${attempt}/3); retrying in $((attempt * 10))s" >&2
    sleep $((attempt * 10))
  done
  echo "could not reach ${url} after 3 attempts" >&2
  return 1
}

if [[ "$VERSION" == "latest-release" ]]; then
  echo "resolving latest Exacto release..."
  release_json="$(gh_api "https://api.github.com/repos/${api_repo}/releases/latest")"
  # strict=False: release notes routinely contain raw control characters, which
  # the default JSON decoder rejects.
  read -r tag asset_url <<<"$(printf '%s' "$release_json" | python3 -c '
import json, sys
release = json.loads(sys.stdin.read(), strict=False)
tarballs = [a for a in release["assets"] if a["name"].endswith(".tar.gz")]
url = tarballs[0]["browser_download_url"] if tarballs else release["tarball_url"]
print(release["tag_name"], url)')"

  echo "installing Exacto ${tag} from ${asset_url}"
  curl -fsSL "$asset_url" -o "$WORK/exacto.tar.gz"
  python -m pip install "$WORK/exacto.tar.gz" --verbose
  source_ref="$tag"
  source_kind="release"
else
  echo "building Exacto from ref ${VERSION}"
  rm -rf "$WORK/checkout"
  git clone --depth 1 --branch "$VERSION" "$REPO" "$WORK/checkout"
  python -m pip install "$WORK/checkout" --verbose
  source_ref="$(git -C "$WORK/checkout" rev-parse --short HEAD)"
  source_kind="source"
fi

installed="$(python -c 'import importlib.metadata as m; print(m.version("exacto"))')"
echo "exacto ${installed} installed"
exacto --help > /dev/null

python - "$installed" "$source_ref" "$source_kind" "$RESULTS/environment.json" <<'PY'
import json
import platform
import subprocess
import sys

installed, source_ref, source_kind, out_path = sys.argv[1:5]


def version_of(command):
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        return (result.stdout or result.stderr).strip().splitlines()[0]
    except FileNotFoundError:
        return None


payload = {
    "exacto_version": installed,
    "exacto_ref": source_ref,
    "exacto_source": source_kind,
    "python": sys.version.split()[0],
    "platform": platform.platform(),
    "samtools": version_of(["samtools", "--version"]),
    "minimap2": version_of(["minimap2", "--version"]),
    "rnabloom": version_of(["rnabloom", "-v"]),
    # Every tool a method in pipeline/methods.py can invoke. Recorded whether
    # or not this leg uses it, so the site reports the environment the run had
    # rather than the subset one arm happened to touch. Absent tools record as
    # null rather than failing: a leg that does not need isONform should not
    # die because it is missing.
    "rnaspades": version_of(["rnaspades.py", "--version"]),
    "isonclust": version_of(["isONclust", "--version"]),
    "isoncorrect": version_of(["run_isoncorrect", "--version"]),
    "isonform": version_of(["isONform_parallel", "--version"]),
}

with open(out_path, "w") as handle:
    json.dump(payload, handle, indent=2)
    handle.write("\n")

print(json.dumps(payload, indent=2))
PY
