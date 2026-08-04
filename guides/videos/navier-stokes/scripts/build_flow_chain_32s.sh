#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
video_dir="$(cd -- "$script_dir/.." && pwd)"
chain_dir="$video_dir/renders/flow/chain"
output="$video_dir/renders/edit/navier-stokes-flow-chain-v01-32s.mp4"

s01="$video_dir/renders/flow/drafts/S01_flow_lite_720p_take01.mp4"
c01="$chain_dir/C01_flow_lite_first-last_8s_take01.mp4"
c02="$chain_dir/C02_flow_lite_first-last_8s_take01.mp4"
c03="$chain_dir/C03_flow_lite_first-last_8s_take01.mp4"

for input in "$s01" "$c01" "$c02" "$c03"; do
  if [[ ! -f "$input" ]]; then
    echo "Missing required input: $input" >&2
    exit 1
  fi
done

mkdir -p "$(dirname -- "$output")"

ffmpeg -hide_banner -y \
  -i "$s01" -i "$c01" -i "$c02" -i "$c03" \
  -filter_complex \
  "[0:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v0];
   [1:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v1];
   [2:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v2];
   [3:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v3];
   [0:a]aformat=sample_rates=48000:channel_layouts=stereo,asetpts=PTS-STARTPTS[a0];
   [1:a]aformat=sample_rates=48000:channel_layouts=stereo,asetpts=PTS-STARTPTS[a1];
   [2:a]aformat=sample_rates=48000:channel_layouts=stereo,asetpts=PTS-STARTPTS[a2];
   [3:a]aformat=sample_rates=48000:channel_layouts=stereo,asetpts=PTS-STARTPTS[a3];
   [v0][a0][v1][a1][v2][a2][v3][a3]concat=n=4:v=1:a=1[v][a]" \
  -map "[v]" -map "[a]" \
  -c:v libx264 -preset medium -crf 18 -pix_fmt yuv420p \
  -c:a aac -b:a 192k -ar 48000 \
  -movflags +faststart \
  "$output"

echo "Created: $output"
