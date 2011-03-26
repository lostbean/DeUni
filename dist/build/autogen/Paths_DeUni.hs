module Paths_DeUni (
    version,
    getBinDir, getLibDir, getDataDir, getLibexecDir,
    getDataFileName
  ) where

import Data.Version (Version(..))
import System.Environment (getEnv)

version :: Version
version = Version {versionBranch = [0,0,1], versionTags = []}

bindir, libdir, datadir, libexecdir :: FilePath

bindir     = "/home/edgar/.cabal/bin"
libdir     = "/home/edgar/.cabal/lib/DeUni-0.0.1/ghc-7.0.2"
datadir    = "/home/edgar/.cabal/share/DeUni-0.0.1"
libexecdir = "/home/edgar/.cabal/libexec"

getBinDir, getLibDir, getDataDir, getLibexecDir :: IO FilePath
getBinDir = catch (getEnv "DeUni_bindir") (\_ -> return bindir)
getLibDir = catch (getEnv "DeUni_libdir") (\_ -> return libdir)
getDataDir = catch (getEnv "DeUni_datadir") (\_ -> return datadir)
getLibexecDir = catch (getEnv "DeUni_libexecdir") (\_ -> return libexecdir)

getDataFileName :: FilePath -> IO FilePath
getDataFileName name = do
  dir <- getDataDir
  return (dir ++ "/" ++ name)
