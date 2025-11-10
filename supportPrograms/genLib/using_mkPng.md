This is not complete. Will finish when I get around to
  getting mkPng working.

# mkPng:

Not yet working. When working, it will not be very good
  since it does not apply compression to output PNGs.
  When working code to make a PNG.

Uses st\_mkPng structure. You should not fiddle with this
  structure.

Inorder to get a PNG structure reader for use you need
  to blank the structure with `blank_st_mkPng`. Then
  initialize with `init_st_mkPng`. Next you allocate
  memory and set the image dimentions
  with `setup_st_mkPng`.

You can skip the initial steps by calling `mk_st_png`.

To reuse a png structure, call `blank_st_mkPng` and then
  `setup_st_mkPng`.

You can add or remove colors to your pallete
  with `addCol_st_mkPng` and `rmCol_st_mkPng`.

You can graph a pixel with `addPixel_st_mkPng` or a
  straight line with `addBar_st_mkPng`.

When finished print your png with `print_st_mkPng`. Then
  free you png struct with `freeStack_st_mkPng` (if on
  stack) or `freeHeap_st_mkPng` if on heap.

