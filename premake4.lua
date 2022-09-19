
local action = _ACTION or ""

solution "libtess2"
	location ( "Build" )
	configurations { "Debug", "Release" }
	platforms {"native", "x64", "x32"}
  
	configuration "Debug"
		defines { "DEBUG" }
		flags { "Symbols", "ExtraWarnings"}

	configuration "Release"
		defines { "NDEBUG" }
		flags { "Optimize", "ExtraWarnings"}    


	project "tess2"
		language "C"
		kind "StaticLib"
		includedirs { "Include", "Source" }
		files { "Source/*.c" }
		targetdir("Build")

	-- more dynamic example
	project "example"
		kind "ConsoleApp"
		language "C"
		links { "tess2" }
		files { "Example/example.c", "Contrib/*.c" }
		includedirs { "Include", "Contrib" }
		targetdir("Build")
	 
		configuration { "linux" }
			 linkoptions { "`pkg-config --libs glfw3`" }
			 links { "GL", "GLU", "m", "GLEW" }
			 defines { "NANOVG_GLEW" }

		configuration { "windows" }
			 links { "glfw3", "gdi32", "winmm", "user32", "GLEW", "glu32","opengl32" }
			 defines { "NANOVG_GLEW" }

		configuration { "macosx" }
			links { "glfw3" }
			linkoptions { "-framework OpenGL", "-framework Cocoa", "-framework IOKit", "-framework CoreVideo" }


	-- more dynamic example
	project "emboss"
		kind "ConsoleApp"
		language "C"
		links { "tess2" }
		files { "Emboss/emboss.c" }
		includedirs { "Include", "Contrib" }
		targetdir("Build")

		configuration { "macosx" }
			-- links { "glfw3" }
			-- linkoptions { "-framework OpenGL", "-framework Cocoa", "-framework IOKit", "-framework CoreVideo" }


	-- more dynamic example
	project "wemboss"
		kind "ConsoleApp"
		language "C"
		links { "tess2" }
		files { "WEmboss/main.c", "Contrib/frozen.c", "WEmboss/triangles.c", "WEmboss/surface.c" }
		includedirs { "Include", "Contrib", "WEmboss", "Source", "/opt/homebrew/Cellar/glfw/3.3.7/include" }
		targetdir("Build")

		configuration { "macosx" }
			-- links { "glfw3" }
			-- linkoptions { "-framework OpenGL", "-framework Cocoa", "-framework IOKit", "-framework CoreVideo" }


	-- more dynamic example
	project "geowrap"
	kind "ConsoleApp"
	language "C"
	links { "tess2" }
	files { "geowrap/main.c", "Contrib/frozen.c", "WEmboss/triangles.c", "WEmboss/surface.c" }
	includedirs { "Include", "Contrib", "geowrap", "Source", "WEmboss", "/opt/homebrew/include" }
	targetdir("Build")

	configuration { "macosx" }
		-- links { "glfw3" }
		-- linkoptions { "-framework OpenGL", "-framework Cocoa", "-framework IOKit", "-framework CoreVideo" }


	-- more dynamic example
	project "jtest"
		kind "ConsoleApp"
		language "C"
		links { "tess2" }
		files { "jtest/main.c", "Contrib/frozen.c", "WEmboss/surface.c" }
		includedirs { "Include", "Contrib", "WEmboss", "Source", "/opt/homebrew/Cellar/glfw/3.3.7/include"}
		targetdir("Build")

		configuration { "macosx" }
			-- links { "glfw3" }
			-- linkoptions { "-framework OpenGL", "-framework Cocoa", "-framework IOKit", "-framework CoreVideo" }

			
