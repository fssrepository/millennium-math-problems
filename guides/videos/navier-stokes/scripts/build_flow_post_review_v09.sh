#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
guide_dir="$(cd -- "${script_dir}/.." && pwd)"

# Start from V07, before the incorrectly placed V08 gain correction.
default_input="${guide_dir}/renders/edit/navier-stokes-flow-post-v07-review-720p.mp4"
input_file="${1:-${default_input}}"
output_file="${2:-${guide_dir}/renders/edit/navier-stokes-flow-post-v09-review-720p.mp4}"

if [[ ! -f "${input_file}" ]]; then
  echo "Missing required file: ${input_file}" >&2
  exit 1
fi

mkdir -p -- "$(dirname -- "${output_file}")"

# The recurring pulse that builds up audibly is the 7.67 s occurrence, not
# the similar one at 6.67 s. Preserve the video bit-for-bit and leave V07's
# two-millisecond seam repair unchanged. Apply only a smooth 120 ms gain ride
# at 7.625-7.745 s, reaching -4 dB around 7.67-7.69 s.
ffmpeg -hide_banner -y \
  -i "${input_file}" \
  -map 0:v:0 -map 0:a:0 \
  -c:v copy \
  -af "volume=eval=frame:volume='if(lt(t,7.625),1,if(lt(t,7.665),1-(1-0.630957)*(t-7.625)/0.040,if(lt(t,7.695),0.630957,if(lt(t,7.745),0.630957+(1-0.630957)*(t-7.695)/0.050,1))))'" \
  -c:a aac -b:a 256k -ar 48000 \
  -t 29.700 -movflags +faststart \
  "${output_file}"

ffprobe -v error \
  -show_entries format=duration,size:stream=index,codec_name,width,height,r_frame_rate,sample_rate,channels \
  -of default=noprint_wrappers=1 \
  "${output_file}"
