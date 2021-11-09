#!/bin/env fish
g++ ./src/**.cc (cat ./buildOptions.fish) -o ./bin/a.out
./bin/a.out
