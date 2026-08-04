#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
video_dir="$(cd -- "$script_dir/.." && pwd)"
assets="$video_dir/assets"
output="$video_dir/renders/edit/navier-stokes-five-anchor-morph-previs-v01-11s.mp4"

k0="$assets/flow-chain/chain-anchor-start-s01-first.png"
k1="$assets/flow-chain/chain-anchor-00-s01-last.png"
k2="$assets/flow-chain/chain-anchor-01-s02-last.png"
k3="$assets/flow-chain/chain-anchor-02-s05-last-clean.png"
k4="$assets/keyframe-06-frontier-heightfield-v2.png"

for input in "$k0" "$k1" "$k2" "$k3" "$k4"; do
  if [[ ! -f "$input" ]]; then
    echo "Missing required keyframe: $input" >&2
    exit 1
  fi
done

mkdir -p "$(dirname -- "$output")"

ffmpeg -hide_banner -y \
  -loop 1 -t 3 -i "$k0" \
  -loop 1 -t 3 -i "$k1" \
  -loop 1 -t 3 -i "$k2" \
  -loop 1 -t 3 -i "$k3" \
  -loop 1 -t 3 -i "$k4" \
  -filter_complex \
  "[0:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[v0];
   [1:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[v1];
   [2:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[v2];
   [3:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[v3];
   [4:v]scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=24,format=yuv420p,setpts=PTS-STARTPTS[v4];
   [v0][v1]xfade=transition=fade:duration=1:offset=2[x1];
   [x1][v2]xfade=transition=fade:duration=1:offset=4[x2];
   [x2][v3]xfade=transition=fade:duration=1:offset=6[x3];
   [x3][v4]xfade=transition=fade:duration=1:offset=8[v]" \
  -map "[v]" -t 11 \
  -c:v libx264 -preset medium -crf 18 -pix_fmt yuv420p \
  -movflags +faststart \
  "$output"

echo "Created: $output"
