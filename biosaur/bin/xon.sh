#!/bin/sh

# set locale if it is totally undefined
if [ -z "${LC_ALL+x}" ] && [ -z "${LC_CTYPE+x}" ] && \
   [ -z "${LANG+x}" ] && [ -z "${LANGUAGE+x}" ]; then
  export LANG=C.UTF-8
fi

# run python
exec "/home/zloydanny/working_area_latest/bio/bio_project/biosaur/biosaur/bin/python3" -u -m xonsh "$@"
