#!/bin/sh
# The checks every change must pass. Each backend has its own
# elementary-function library, so a tolerance that holds on one target can
# fail on another; all four are tested.
set -eu
cd "$(dirname "$0")/.."
moon fmt --check
# The documentation rules shared by every gaato MoonBit repository; the
# self-test first checks that the audit still catches what it should.
moon run scripts/docs_audit.mbtx --self-test
moon check --target all
for target in wasm-gc js native wasm; do
  moon test --target "$target"
done
before=$(find src -name pkg.generated.mbti -exec sha256sum {} + | sort)
moon info
after=$(find src -name pkg.generated.mbti -exec sha256sum {} + | sort)
if [ "$before" != "$after" ]; then
  echo "pkg.generated.mbti changed: the public API moved; review the change" >&2
  exit 1
fi
