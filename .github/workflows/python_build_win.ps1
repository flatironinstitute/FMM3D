$ErrorActionPreference = "Stop"
Set-Variable -Name PYTHON -Value (Get-Command python).definition
Set-Variable -Name MSYSTEM -Value MINGW64

# setup setup.cfg
New-Item -Force -Path .\python -Name "setup.cfg" -ItemType "file" -Value "[build]`r`ncompiler=mingw32`r`n[build_ext]`r`ncompiler=mingw32"

# Setup the make.inc file
Copy-Item -Path make.inc.windows.mingw -Destination make.inc
Add-Content -Path make.inc -Value "PYTHON=""$PYTHON"""
Add-Content -Path make.inc -Value "FFLAGS+= -fallow-argument-mismatch -march=x86-64"
Add-Content -Path make.inc -Value "CFLAGS+= -march=x86-64"
Add-Content -Path make.inc -Value "CXXFLAGS+= -march=x86-64"

# Setup the distutils.cfg file
Set-Variable distutils_cfg -Value ([IO.Path]::Combine((Split-Path -Path $PYTHON), "Lib", 'distutils', 'distutils.cfg'))
Set-Content -Path $distutils_cfg -Value "[build]`r`ncompiler=mingw32`r`n[build_ext]`r`ncompiler=mingw32"
python -m pip install --upgrade setuptools wheel numpy pip
if (-not $?) {throw "Failed pip install"}

# call make
Set-Variable repo_root -Value ([IO.Path]::Combine($PSScriptRoot, '..', '..'))
c:\msys64\usr\bin\env MSYSTEM=MINGW64 c:\msys64\usr\bin\bash.exe -lc "cd '$repo_root' && make python-dist"
if (-not $?) {throw "Failed make python-dist"}

# Move the required DLLs inside the wheel
wheel.exe unpack (get-item .\python\dist\*.whl).FullName -d .\tpm
if (-not $?) {throw "Failed unpack wheel"}
Set-Variable unpacked_wheel -Value (get-item .\tpm\fmm3dpy-*).FullName
Copy-Item -Path C:\msys64\mingw64\bin\libgcc_s_seh-*.dll -Destination ([IO.Path]::Combine($unpacked_wheel, 'fmm3dpy'))
Copy-Item -Path C:\msys64\mingw64\bin\libgomp-*.dll -Destination ([IO.Path]::Combine($unpacked_wheel, 'fmm3dpy'))
Copy-Item -Path C:\msys64\mingw64\bin\libgfortran-*.dll -Destination ([IO.Path]::Combine($unpacked_wheel, 'fmm3dpy'))
Copy-Item -Path C:\msys64\mingw64\bin\libwinpthread-*.dll -Destination ([IO.Path]::Combine($unpacked_wheel, 'fmm3dpy'))
Copy-Item -Path C:\msys64\mingw64\bin\libquadmath-*.dll -Destination ([IO.Path]::Combine($unpacked_wheel, 'fmm3dpy'))
New-Item -Path .\wheelhouse -ItemType Directory -Force
wheel.exe pack $unpacked_wheel -d .\wheelhouse
if (-not $?) {throw "Failed pack wheel"}

# Cleanup
Remove-Item -Path $unpacked_wheel -Force -Recurse
Remove-Item -Path .\python\build -Force -Recurse
Remove-Item -Path .\python\dist -Force -Recurse
