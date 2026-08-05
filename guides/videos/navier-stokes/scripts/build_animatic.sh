#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"
assets_dir="${guide_dir}/assets"
drafts_dir="${guide_dir}/renders/flow/drafts"
edits_dir="${guide_dir}/renders/flow/edits"
output_dir="${guide_dir}/renders/edit"
overlay_file="${guide_dir}/04_editorial/overlays/animatic-overlays.ass"
output_file="${1:-${output_dir}/navier-stokes-animatic-v02-raw-plates-720x1280.mp4}"
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

require_file "${assets_dir}/keyframe-01-world.png"
require_file "${assets_dir}/keyframe-03-checkpoint.png"
require_file "${assets_dir}/keyframe-06-frontier.png"
require_file "${drafts_dir}/S02_flow_lite_720p_take02.mp4"
require_file "${drafts_dir}/S05_flow_lite_720p_take01.mp4"
require_file "${overlay_file}"

render_still() {
  local input_file="$1"
  local frame_count="$2"
  local output_path="$3"

  ffmpeg -loglevel error -y \
    -framerate 30 -loop 1 -i "${input_file}" \
    -vf "scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=30,trim=end_frame=${frame_count},setpts=PTS-STARTPTS,format=yuv420p" \
    -frames:v "${frame_count}" \
    -an -c:v libx264 -preset veryfast -crf 18 \
    "${output_path}"
}

render_video() {
  local input_file="$1"
  local frame_count="$2"
  local output_path="$3"

  ffmpeg -loglevel error -y \
    -i "${input_file}" \
    -vf "scale=720:1280:force_original_aspect_ratio=increase,crop=720:1280,fps=30,trim=end_frame=${frame_count},setpts=PTS-STARTPTS,format=yuv420p" \
    -frames:v "${frame_count}" \
    -an -c:v libx264 -preset veryfast -crf 18 \
    "${output_path}"
}

# Exact 30 fps allocation: 105+240+27+144+115+125+135 = 891 frames = 29.70 s.
# The full raw eight-second S02 source spans all of S02 and the first two
# seconds of S03, preserving the generated ending instead of freezing it.
# S05 is also shown at normal source speed; no salvage retiming is used here.
render_still "${assets_dir}/keyframe-01-world.png" 105 "${work_dir}/s01.mp4"
render_video "${drafts_dir}/S02_flow_lite_720p_take02.mp4" 240 "${work_dir}/s02-s03a.mp4"
render_still "${assets_dir}/keyframe-03-checkpoint.png" 27 "${work_dir}/s03b.mp4"
render_still "${assets_dir}/keyframe-06-frontier.png" 144 "${work_dir}/s04.mp4"
render_video "${drafts_dir}/S05_flow_lite_720p_take01.mp4" 115 "${work_dir}/s05.mp4"
render_still "${assets_dir}/keyframe-06-frontier.png" 125 "${work_dir}/s06.mp4"
render_still "${assets_dir}/keyframe-06-frontier.png" 135 "${work_dir}/s07.mp4"

ffmpeg -loglevel error -y \
  -i "${work_dir}/s01.mp4" \
  -i "${work_dir}/s02-s03a.mp4" \
  -i "${work_dir}/s03b.mp4" \
  -i "${work_dir}/s04.mp4" \
  -i "${work_dir}/s05.mp4" \
  -i "${work_dir}/s06.mp4" \
  -i "${work_dir}/s07.mp4" \
  -f lavfi -t 29.70 -i "anullsrc=r=48000:cl=stereo" \
  -filter_complex "[0:v][1:v][2:v][3:v][4:v][5:v][6:v]concat=n=7:v=1:a=0[sequence];[sequence]drawbox=x=0:y=0:w=iw:h=ih:color=black@0.58:t=fill:enable='between(t,7.4,9.5)+gte(t,25.2)',subtitles='${overlay_file}'[video]" \
  -map "[video]" -map 7:a \
  -frames:v 891 -t 29.70 \
  -c:v libx264 -preset medium -crf 18 -profile:v high -level 4.0 \
  -c:a aac -b:a 128k -ar 48000 \
  -pix_fmt yuv420p -movflags +faststart \
  "${output_file}"

ffprobe -v error \
  -show_entries format=duration,size:stream=index,codec_type,codec_name,width,height,pix_fmt,r_frame_rate,sample_rate,channels \
  -of default=noprint_wrappers=1 \
  "${output_file}"
