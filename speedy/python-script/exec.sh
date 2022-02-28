#!/bin/bash
sudo docker run -it -e DISPLAY=$DISPLAY --net host -v /tmp/.X11-unix:/tmp/.X11-unix -v $HOME/.Xauthority:/root/.Xauthority -v /data02/oishi/:/home/ --shm-size 16g oishi/cdo
# /bin/bash -c "cd /home/speedy-epyc/speedy/python-script; python3 rmse.py convert"
