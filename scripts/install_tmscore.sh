
# Update the package lists
sudo apt-get update

# Install g++
sudo apt-get install -y g++

# Download TMscore.cpp
wget https://seq2fun.dcmb.med.umich.edu/TM-score/TMscore.cpp -O TMscore.cpp

# Check if the download was successful
if [ $? -ne 0 ]; then
    echo "Failed to download TMscore.cpp"
    exit 1
fi

# Compile the TMscore.cpp with the specified flags
g++ -static -O3 -ffast-math -lm -o TMscore TMscore.cpp

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful. The executable is named TMscore."
else
    echo "Compilation failed."
fi