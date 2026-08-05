#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"

default_input="/home/raxim/Downloads/Initial_Scene_-_2026-08-05_202608051411.mp4"
input_file="${1:-${default_input}}"
output_file="${2:-${guide_dir}/renders/edit/navier-stokes-flow-post-v01-29.70s.mp4}"
overlay_file="${guide_dir}/04_editorial/overlays/flow-final-overlays-v01.ass"

for required_file in "${input_file}" "${overlay_file}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "Missing required file: ${required_file}" >&2
    exit 1
  fi
done

mkdir -p -- "$(dirname -- "${output_file}")"

# The downloaded Flow chain is 31.000 seconds because each Extend reuses about
# one second of context. The visible discontinuity begins at source frame 552
# (23.000 s), and the new geometry then persists. Drop only source frames 552
# and 553, bridge the adjacent states with a minimal two-frame dissolve, and
# use the stable generated tail to retain the locked 29.700-second runtime.
# Audio follows the identical edit, then receives transparent impulse cleanup.
ffmpeg -hide_banner -y \
  -i "${input_file}" \
  -filter_complex "\
[0:v]split=2[va0][vb0];\
[va0]trim=start=0:end=23.000,setpts=PTS-STARTPTS[va];\
[vb0]trim=start=23.083333:end=29.866667,setpts=PTS-STARTPTS[vb];\
[va][vb]xfade=transition=fade:duration=0.083333:offset=22.916667,\
fps=30,scale=1080:1920:flags=lanczos,setsar=1,\
tpad=stop_mode=clone:stop_duration=0.034,trim=duration=29.700,\
subtitles='${overlay_file}',format=yuv420p[vout];\
[0:a]aresample=48000,asplit=6[a0][a1][a2][a3][a4][a5];\
[a0]atrim=start=0:end=6.880,asetpts=PTS-STARTPTS[p0];\
[a1]atrim=start=6.880:end=7.120,asetpts=PTS-STARTPTS,adeclick[p1];\
[a2]atrim=start=7.120:end=14.880,asetpts=PTS-STARTPTS[p2];\
[a3]atrim=start=14.880:end=15.120,asetpts=PTS-STARTPTS,adeclick[p3];\
[a4]atrim=start=15.120:end=23.000,asetpts=PTS-STARTPTS[aa];\
[a5]atrim=start=23.083333:end=29.866667,asetpts=PTS-STARTPTS[ab];\
[aa][ab]acrossfade=d=0.083333:c1=tri:c2=tri[tail];\
[p0][p1][p2][p3][tail]concat=n=5:v=0:a=1,\
apad=pad_dur=0.050,atrim=duration=29.700,volume=-1.2dB[aout]" \
  -map "[vout]" -map "[aout]" \
  -frames:v 891 -t 29.700 \
  -c:v libx264 -preset slow -crf 17 -profile:v high -level 4.2 \
  -c:a aac -b:a 320k -ar 48000 \
  -movflags +faststart \
  "${output_file}"

ffprobe -v error \
  -show_entries format=duration,size:stream=index,codec_name,width,height,r_frame_rate,sample_rate,channels \
  -of default=noprint_wrappers=1 \
  "${output_file}"
