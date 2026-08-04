#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"
input_file="${1:-${guide_dir}/assets/keyframe-06-frontier.png}"
output_file="${2:-${guide_dir}/renders/flow/edits/S06_keyframe_4p16s_dolly_silent.mp4}"

if [[ ! -f "${input_file}" ]]; then
  echo "Missing input: ${input_file}" >&2
  exit 1
fi

mkdir -p -- "$(dirname -- "${output_file}")"

# The Flow take generated visible 4%/6% labels and morphed checkpoint cubes
# into spike-like glyphs. Animate the audited canonical frame instead: begin
# slightly closer, ease out to the exact full reference, and omit source audio.
ffmpeg -loglevel error -y \
  -framerate 30 \
  -loop 1 \
  -i "${input_file}" \
  -filter_complex "[0:v]scale=720:1280,split=2[movein][holdin];[movein]zoompan=z='1.04-0.04*on/111':x='iw/2-(iw/zoom/2)':y='ih/2-(ih/zoom/2)':d=112:s=720x1280:fps=30,trim=end_frame=112,setpts=PTS-STARTPTS[move];[holdin]fps=30,trim=end_frame=13,setpts=PTS-STARTPTS[hold];[move][hold]concat=n=2:v=1:a=0,format=yuv420p[outv]" \
  -map "[outv]" \
  -c:v libx264 \
  -preset medium \
  -crf 18 \
  -an \
  -movflags +faststart \
  "${output_file}"

ffprobe -v error \
  -show_entries format=duration,size:stream=index,codec_type,codec_name,width,height,pix_fmt,r_frame_rate \
  -of default=noprint_wrappers=1 \
  "${output_file}"
