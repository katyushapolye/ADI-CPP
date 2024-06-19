N=16
FILELOGGING=1 #se vai exportar quadros e log de error
FRAMES=5 #framerate do video final -> buscar algo em torno de 30quadros por segundo pra ficar legal
NAME=graph.mp4
rm Frames/* 
rm Fields/*
rm -r bin/solver
rm VectorFields/*
rm VectorFrames/*
rm $NAME

g++ -m64 -mpc64 -mfpmath=sse -std=c++14 -Wall -Wno-sign-compare -o bin/solver src/source/* -I${workspaceFolder}/src/headers


./bin/solver $N $FILELOGGING

#python3 bin/plotter.py $N 

python3 bin/plotter_vec.py $N 

ffmpeg -hide_banner -loglevel error -framerate $FRAMES -i 'VectorFrames/VectorFrame_%d.png' -c:v libx264 -pix_fmt yuv420p $NAME

mpv $NAME
