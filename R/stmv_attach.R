
stmv_attach = function( method, pointer ) {
  #\\ generic method to attach data pointer from bigmemory or ff storage pointers
  return(
    switch( method,
      bigmemory.ram = bigmemory::attach.big.matrix(pointer),
      bigmemory.filebacked = bigmemory::attach.big.matrix(pointer),
      ff=pointer )
  )
}

