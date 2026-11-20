/*
  CRAM move prior for Stockfish (experimental add-on, not part of upstream Stockfish)

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
*/

// Standalone CRAM net embedding for universal binary builds. The net is only about 100 KB, so
// unlike the NNUE net it is simply embedded into every slice, including the macOS x86 one.
// See cram.cpp for how the symbols are used.

extern const unsigned char gEmbeddedCRAMData[] =
#ifdef __has_embed
{
    #embed "model-i8.bin"
};
const unsigned int cramPadding = 0;
#else
    #include "cram_dump.inc"
  ;
const unsigned int cramPadding = 1;  // trailing NUL byte
#endif
extern const unsigned int gEmbeddedCRAMSize = sizeof(gEmbeddedCRAMData) - cramPadding;
