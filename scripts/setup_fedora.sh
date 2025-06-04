#!/bin/bash

dnf install -y powerline powerline-fonts cmake clang build-ninja fmt-devel qt6-qtbase-devel cairo-devel tbb-devel

cat << EOL >> ~/.bashrc
export CMAKE_GENERATOR=Ninja
if [ -f \$(which powerline-daemon) ]; then
  powerline-daemon -q
  POWERLINE_BASH_CONTINUATION=1
  POWERLINE_BASH_SELECT=1
  . /usr/share/powerline/bash/powerline.sh
fi
EOL
