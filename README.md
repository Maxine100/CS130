# CS130
Here's some code for CS130, which I took in the winter of 2020.

## Ray Tracer
Use scons to make a make file.
```
$ scons
```
To create an image, use
```
$ ./ray_tracer -i <##>.txt
```
That will generate an image into output.png.

To comopare, use
```
$ ./ray_tracer -i <##>.txt -s <##>.png
```
That will compare output.png and the image that's supposed to be generated. The differences will be in diff.png.

To run the grading script, use
```
$ ./grading-script.py .
```

## OpenGL Graphics Pipeline
To compile this project, use the same commands as those for the ray tracer, but replace ray_tracer with driver.
