set -e

## compile
javac -d bin -sourcepath src -cp lib/jscience.jar src/Test.java

## run
java -cp .:bin:lib/jscience.jar Test
