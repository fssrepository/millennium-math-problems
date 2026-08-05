#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"
drafts_dir="${guide_dir}/renders/flow/drafts"
edits_dir="${guide_dir}/renders/flow/edits"
assets_dir="${guide_dir}/assets"
output_dir="${guide_dir}/renders/edit"
overlay_file="${guide_dir}/04_editorial/overlays/raw-flow-review.ass"
output_file="${1:-${output_dir}/navier-stokes-raw-flow-review-v01-44p2s.mp4}"
work_dir="$(mktemp -d)"

cleanup() {
  rm -rf -- "${work_dir}"
}
trap cleanup EXIT

mkdir -p -- "${output_dir}" "$(dirname -- "${output_file}")"

require_file() {
  if [[ ! -f "$1" ]]; then
    echo "Missing required input: $1" >&2
    exit 1
  fi
}

s01="${drafts_dir}/S01_flow_lite_720p_take01.mp4"
s02="${drafts_dir}/S02_flow_lite_720p_take02.mp4"
s05="${drafts_dir}/S05_flow_lite_720p_take01.mp4"
s06="${edits_dir}/S06_flow_lite_take01_8s_label_removed.mp4"

require_file "${s01}"
require_file "${s02}"
require_file "${s05}"
require_file "${s06}"
require_file "${assets_dir}/keyframe-03-checkpoint.png"
require_file "${assets_dir}/keyframe-06-frontier.png"
require_file "${overlay_file}"

render_still_with_silence() {
  local input_file="$1"
  local duration="$2"
  local output_path="$3"

  ffmpeg -loglevel error -y \
    -framerate 24 -loop 1 -i "${input_file}" \
    -f lavfi -i "anullsrc=r=48000:cl=stereo" \
    -t "${duration}" \
    -vf "scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,format=yuv420p" \
    -c:v libx264 -preset veryfast -crf 18 -r 24 \
    -c:a aac -b:a 192k -ar 48000 -shortest \
    "${output_path}"
}

render_still_with_silence "${assets_dir}/keyframe-03-checkpoint.png" 2.9 "${work_dir}/s03.mp4"
render_still_with_silence "${assets_dir}/keyframe-06-frontier.png" 4.8 "${work_dir}/s04.mp4"
render_still_with_silence "${assets_dir}/keyframe-06-frontier.png" 4.5 "${work_dir}/s07.mp4"

ffmpeg -loglevel error -y \
  -i "${s01}" \
  -i "${s02}" \
  -i "${work_dir}/s03.mp4" \
  -i "${work_dir}/s04.mp4" \
  -i "${s05}" \
  -i "${s06}" \
  -i "${work_dir}/s07.mp4" \
  -filter_complex "[0:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v0];[0:a]aresample=48000,asetpts=PTS-STARTPTS[a0];[1:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v1];[1:a]aresample=48000,asetpts=PTS-STARTPTS[a1];[2:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v2];[2:a]aresample=48000,asetpts=PTS-STARTPTS[a2];[3:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v3];[3:a]aresample=48000,asetpts=PTS-STARTPTS[a3];[4:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v4];[4:a]aresample=48000,asetpts=PTS-STARTPTS[a4];[5:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v5];[5:a]aresample=48000,asetpts=PTS-STARTPTS[a5];[6:v]fps=24,scale=720:1280,setsar=1,setpts=PTS-STARTPTS[v6];[6:a]aresample=48000,asetpts=PTS-STARTPTS[a6];[v0][a0][v1][a1][v2][a2][v3][a3][v4][a4][v5][a5][v6][a6]concat=n=7:v=1:a=1[basev][audio];[basev]subtitles='${overlay_file}'[video]" \
  -map "[video]" -map "[audio]" \
  -t 44.2 \
  -c:v libx264 -preset medium -crf 18 -profile:v high -level 4.0 \
  -c:a aac -b:a 192k -ar 48000 \
  -pix_fmt yuv420p -movflags +faststart \
  "${output_file}"

ffprobe -v error \
  -show_entries format=duration,size:stream=index,codec_type,codec_name,width,height,pix_fmt,r_frame_rate,sample_rate,channels \
  -of default=noprint_wrappers=1 \
  "${output_file}"
