#!/bin/bash
if [ -f "wavefield/bwd_merged/merged_instaseis_db.nc4" ]; then
  echo "Wavefields available from cache"
else
  echo "Wavefields need to be downloaded"

  # Check, whether encrypted token for Google Drive account is available
  # Local gdrive installation needs to be activated with it, but since it
  # offers read/write access, it should not be shared with other accounts 
  # than github/seismology.
  # Download from gdrive is much faster (50-100MB/s) than from LMU (<10 MB/s),
  # but pull requests from other github accounts than seismology will not have
  # access to it. 
  wget -q https://polybox.ethz.ch/index.php/s/35W6wls9ARpUHiD/download -O kerner_wavefields.tar.bz2
  tar -xvf kerner_wavefields.tar.bz2

  # if [ $encrypted_50ebf69bd92e_key == "" ]; then

  #   echo "Downloading from LMU"
  #   # 1. Download archive from LMU server
  #   wget https://www.geophysik.uni-muenchen.de/~staehler/kerner_wavefields.tar.bz2

  #   # 2. Unpack archive
  #   tar -xvf kerner_wavefields.tar.bz2

  # else

  #   echo "Downloading from Google Drive"
  #   # 1. Download gdrive client
  #   wget "https://docs.google.com/uc?id=0B3X9GlR6EmbnQ0FtZmJJUXEyRTA&export=download" -O gdrive
  #   chmod +x ./gdrive

  #   # 2. Decrypt gdrive token
  #   mkdir $HOME/.gdrive
  #   openssl aes-256-cbc -K $encrypted_50ebf69bd92e_key -iv $encrypted_50ebf69bd92e_iv -in TRAVIS/token_v2.json.enc -out $HOME/.gdrive/token_v2.json -d
  #   
  #   # 3. Download wavefield folder with gdrive client
  #   ./gdrive download 0BwU9d5SH6pPQNGpkcW9hWHlpYTA --recursive

  # fi
fi
