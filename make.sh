N=${1:-$DEFAULT_START_N}
FILELOGGING=1 #se vai exportar quadros e log de error
FRAMES=5 #framerate do video final -> buscar algo em torno de 30quadros por segundo pra ficar legal
NAME=graph.mp4
rm Frames/*/* 
rm Fields/*/*
rm -r bin/solver
rm $NAME
rm PressureGradient.mp4

g++ -std=c++14 -Wall -o bin/solver src/source/* -I${workspaceFolder}/src/headers


./bin/solver $N $FILELOGGING

#python3 bin/plotter.py $N 

python3 bin/plotter_vec.py $N 

#plot da aprox
ffmpeg -hide_banner -loglevel error -framerate $FRAMES -i 'Frames/VectorFrames/VectorFrame_%d.png' -c:v libx264 -pix_fmt yuv420p $NAME


#plot da pressao
ffmpeg -hide_banner -loglevel error -framerate $FRAMES -i 'Frames/GradPressureFrames/VectorFrame_%d.png' -c:v libx264 -pix_fmt yuv420p "PressureGradient.mp4"


#plot escalar
#ffmpeg -hide_banner -loglevel error -framerate $FRAMES -i 'Frames/Field_%d.png' -c:v libx264 -pix_fmt yuv420p $NAME

#mpv $NAME
