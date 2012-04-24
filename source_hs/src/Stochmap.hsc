{-# LANGUAGE ForeignFunctionInterface #-}
module Stochmap (calculateAndWrite,calculateStochmap) where

import Foreign
import Foreign.C
import Foreign.Ptr (nullPtr)
import Foreign.C.Types
import Foreign.C.String (CString,newCString)
import System.IO
import System.Posix.Types
import System.Posix
import Data.List
import System.IO.Unsafe

#include "stdio.h"
#include "stochmap.h"

foreign import ccall "" fdopen :: Fd -> CString -> IO (Ptr CFile)
handleToCFile :: Handle -> String -> IO (Ptr CFile)
handleToCFile h m = do iomode <- newCString m
                       fd <- handleToFd h
                       fdopen fd iomode

foreign import ccall "stochmap.h CalculateAndWrite"
     c_CalculateAndWrite :: CInt -> CInt -> CInt -> CInt -> CInt -> Ptr (Ptr (Ptr (Ptr CInt))) -> Ptr (Ptr CInt) -> Ptr CInt -> Ptr CInt -> Ptr (Ptr (Ptr (Ptr (Ptr CDouble)))) -> Ptr (Ptr (Ptr CDouble)) -> Ptr (Ptr CDouble) -> Ptr (Ptr CDouble) -> Ptr CDouble -> Ptr CDouble -> Ptr CFile -> IO (Ptr StochmapResultC)


data StochmapResultC = StochmapResultC {condE :: (Ptr (Ptr (Ptr CDouble))), priorE :: (Ptr (Ptr CDouble))}
instance Storable StochmapResultC where
        sizeOf _ = (#size StochmapResult)
        alignment _ = alignment (undefined::(Ptr CDouble)) -- Ints and Ptrs should both be 4-byte aligned
        peek ptr = do
                condE' <- (#peek StochmapResult, condE) ptr
                priorE' <- (#peek StochmapResult, priorE) ptr
                return StochmapResultC {condE=condE',priorE=priorE'}


parseResult :: Int -> Int -> Int -> StochmapResultC -> IO ([[[Double]]],[[Double]])
parseResult nbranch nproc ncols (StochmapResultC condEC priorEC)= do 
                                                                  priorE'' <- peekArray nbranch priorEC
                                                                  priorE' <- mapM (peekArray nproc) priorE''
                                                                  let priorEConv = map (map realToFrac) priorE'
                                                                  condE''' <- peekArray nbranch condEC
                                                                  condE'' <- mapM (peekArray nproc) condE'''
                                                                  condE' <- (mapM . mapM) (peekArray ncols) condE''
                                                                  let condEConv = (map . map . map) realToFrac condE'
                                                                  return (condEConv,priorEConv)


{-- exposed function --}
--calculateAndWrite :: Int -> Int -> Int -> Int -> Int -> [[Int]] -> [Int] -> [Int] -> [[[[[Double]]]]] -> [[[Double]]] -> [[Double]] -> [[Double]] -> [Double] -> [Double] -> Maybe Handle -> ([[[Double]]],[[Double]])
calculateStochmap nSite nState nBranch nProc nCols lMat multiplicites siteMap partials qset sitelikes pi_i branchLengths mixProbs = System.IO.Unsafe.unsafePerformIO $ calculateAndWrite nSite nState nBranch nProc nCols lMat multiplicites siteMap partials qset sitelikes pi_i branchLengths mixProbs Nothing
calculateAndWrite nSite nState nBranch nProc nCols lMat
            multiplicites sitemap partials qset sitelikes pi_i 
            branchLengths mixProbs handle = do let c_nSite = fromIntegral nSite
                                               let c_nState = fromIntegral nState
                                               let c_nBranch = fromIntegral nBranch
                                               let c_nProc = fromIntegral nProc
                                               let c_nCols = fromIntegral nCols
                                               c_handle <- case handle of
                                                                Just h -> handleToCFile h "w"
                                                                Nothing -> return nullPtr
                                               (c_scales,c_partials) <- makeScalesPtr partials
                                               c_lMat <- newArray2 fromIntegral lMat
                                               c_multiplicites <- newArrayX fromIntegral multiplicites
                                               c_sitemap <- newArrayX fromIntegral sitemap
                                               c_qset <- newArray3 realToFrac qset
                                               c_sitelikes <- newArray2 realToFrac sitelikes
                                               c_pi_i <- newArray2 realToFrac pi_i
                                               c_branchLengths <- newArrayX realToFrac branchLengths
                                               c_mixProbs <- newArrayX realToFrac mixProbs
                                               ans <- c_CalculateAndWrite c_nSite c_nState c_nBranch c_nProc c_nCols c_scales c_lMat c_multiplicites c_sitemap c_partials c_qset c_sitelikes c_pi_i c_branchLengths c_mixProbs c_handle
                                               ans' <- peek ans
                                               parseResult nBranch nProc nSite ans'


newArrayX f xs = newArray $ map f xs

newArray2 :: Storable a => (b->a) -> [[b]] -> IO (Ptr (Ptr a))
newArray2 f xs = newArray =<< (mapM (newArrayX f) xs)
newArray3 f xs = newArray =<< (mapM (newArray2 f) xs)
newArray4 f xs = newArray =<< (mapM (newArray3 f) xs)
newArray5 f xs = newArray =<< (mapM (newArray4 f) xs)


makePartialBranchEnds :: [[[[[Double]]]]] -> IO (Ptr (Ptr (Ptr (Ptr (Ptr CDouble)))))
makePartialBranchEnds = newArray5 realToFrac

makeScales :: [[[[[Double]]]]] -> ([[[[Int]]]],[[[[[Double]]]]])
makeScales xs = umap (umap (umap scale)) xs where
                  umap f x = unzip (map f x)
                  scale m = (scales , m2) where 
                                  scales = map (\y-> head $ reverse $ sort (map makeScale y)) m
                                  m2 = map (\(x,y) -> map (*(10.0**(fromIntegral y))) x ) $ zip m scales

makeScale :: Double -> Int
makeScale d = 0 where -- -1 * (floor $ min 0.0 (5.0 + (log10 d))) where
                  log10 x = (log x) / (log 10)

makeScalesPtr :: [[[[[Double]]]]] -> IO (Ptr (Ptr (Ptr (Ptr CInt))), Ptr (Ptr (Ptr (Ptr (Ptr CDouble)))))
makeScalesPtr xs = do let (scales,mats) = makeScales xs
                      let zipTuple (a,b) = zip a b
                      matsPtr <- makePartialBranchEnds mats
                      scalePtr <- newArray4 fromIntegral scales
                      return (scalePtr,matsPtr)

print1 :: [[[[Int]]]] -> [[[[[Double]]]]] -> IO [[[[()]]]]
print2 :: ([[[Int]]],[[[[Double]]]]) -> IO [[[()]]]
print3 :: ([[Int]],[[[Double]]]) -> IO [[()]]
print4 :: ([Int],[[Double]]) -> IO [()]
print5 :: (Int,[Double]) -> IO ()
print5 (s1,m1) = print $ (show $ length m1) ++ "  Scale " ++ (show s1) ++ " " ++ (show m1)
print4 (s1,m1) = mapM print5 $ zip s1 m1
print3 (s1,m1) = mapM print4 $ zip s1 m1
print2 (s1,m1) = mapM print3 $ zip s1 m1
print1 s m = mapM print2 $ zip s m
