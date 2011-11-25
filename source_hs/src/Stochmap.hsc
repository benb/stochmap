{-# LANGUAGE ForeignFunctionInterface #-}
module Stochmap where

import Foreign
import Foreign.C
import Foreign.C.Types
import Foreign.C.String (CString,newCString)
#include "stochmap.h"

foreign import ccall "stochmap.h CalculateAndWrite"
     c_CalculateAndWrite :: CInt -> CInt -> CInt -> CInt -> CInt -> Ptr (Ptr (Ptr (Ptr CInt))) -> Ptr (Ptr CInt) -> Ptr CInt -> Ptr CInt -> Ptr (Ptr (Ptr (Ptr (Ptr CDouble)))) -> Ptr (Ptr (Ptr CDouble)) -> Ptr (Ptr CDouble) -> Ptr (Ptr CDouble) -> Ptr CDouble -> Ptr CDouble -> Ptr CFile -> IO ()
