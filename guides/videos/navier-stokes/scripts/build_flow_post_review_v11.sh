#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"

default_video_input="/home/raxim/Downloads/Initial_Scene_-_2026-08-05_202608051411.mp4"
default_audio_input="${guide_dir}/renders/edit/navier-stokes-flow-post-v09-review-720p.mp4"
video_input="${1:-${default_video_input}}"
audio_input="${2:-${default_audio_input}}"
output_file="${3:-${guide_dir}/renders/edit/navier-stokes-flow-post-v11-review-720p.mp4}"
overlay_file="${guide_dir}/04_editorial/overlays/flow-final-overlays-v01.ass"

for required_file in "${video_input}" "${audio_input}" "${overlay_file}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "Missing required file: ${required_file}" >&2
    exit 1
  fi
done

mkdir -p -- "$(dirname -- "${output_file}")"
video_tmp="$(mktemp --tmpdir="$(dirname -- "${output_file}")" navier-v11-video-XXXXXX.mp4)"
trap 'rm -f -- "${video_tmp}"' EXIT

# Render the revised visual layer first, then remux the complete approved V09
# audio stream in a separate pass so no audio packets are lost during filtering.
filter_graph="[0:v]split=3[v0][v1][v2];[v0]trim=start=0:end=15.000,setpts=PTS-STARTPTS[va];[v1]trim=start=15.958333:end=23.000,setpts=(PTS-STARTPTS)*1.136095[vb];[v2]trim=start=23.958333:end=31.000,setpts=(PTS-STARTPTS)*1.136095[vc];[va][vb][vc]concat=n=3:v=1:a=0,fps=30,scale=720:1280:flags=lanczos,setsar=1,tpad=stop_mode=clone:stop_duration=0.034,trim=duration=29.700,subtitles='${overlay_file}',format=yuv420p[vout]"

ffmpeg -hide_banner -y -i "${video_input}" -filter_complex "${filter_graph}" -map "[vout]" -frames:v 891 -t 29.700 -an -c:v libx264 -preset medium -crf 18 -profile:v high -level 4.0 -movflags +faststart "${video_tmp}"
ffmpeg -hide_banner -y -i "${video_tmp}" -i "${audio_input}" -map 0:v:0 -map 1:a:0 -c copy -t 29.700 -movflags +faststart "${output_file}"

ffprobe -v error -count_frames -show_entries format=duration,size:stream=index,codec_name,width,height,r_frame_rate,nb_read_frames,start_time,duration,sample_rate,channels -of default=noprint_wrappers=1 "${output_file}"
