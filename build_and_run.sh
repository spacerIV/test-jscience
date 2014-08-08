set -e

## compile
javac -d bin -sourcepath src -cp .:lib/jscience.jar:lib/opencsv-2.4.jar:lib/gson-2.2.4.jar src/Test.java

## run
java -cp .:bin:lib/jscience.jar:lib/opencsv-2.4.jar:lib/gson-2.2.4.jar Test
