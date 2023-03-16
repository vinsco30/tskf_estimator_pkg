#! /bin/bash
awk '/time/{close(f); f="sim" ++c ".txt";} {print > f}' $1
