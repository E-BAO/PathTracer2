# PathTracer2
Continue from PathTracer1, add more complicated materials, environment lights, and depth of field to ray tracer. 

<br/>

## Part 1: Ray Generation and Intersection

### Reflection and refraction

<img src="docs/part1/spheres_256_4_100.png" width="400"/>
With sample rate of 256， 4 samples per light and max depth of 4.

## Part 2: Microfacet Material

<img src="docs/part2/dragon_005_256_4_7.png" width="400"/>
With sample rate of 256， 4 samples per light and max depth of 7.

## Part 3: Environment Light

Importance sampling is introduced to increase the density of sample points in the regions that we are interested in. 

<img src="docs/part3/bunny_env_pdf_4_256_5.png" width="400"/>  <img src="docs/part3/bunny_cu_env_pdf_4_256_5.png" width="400"/>


## Part 4: Depth of Field

Change aperture and focal distance to control the depth of field.

<img src="docs/part5/dragon_256_4_7_06_45.png" width="400"/>  <img src="docs/part5/dragon_256_4_7_04_50.png" width="400"/>
