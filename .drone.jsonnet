local Pipeline(os, arch, version) = {
    kind: "pipeline",
    name: os+" - "+arch+" - Julia "+version,
    platform: {
	os: os,
	arch: arch
    },
    steps: [
	{
	    name: "build",
	    image: "julia:"+version,
	    commands: [
		"julia --project=. --check-bounds=yes --color=yes -e 'using InteractiveUtils; versioninfo(verbose=true); println(split(unsafe_string(ccall(:LLVMGetHostCPUFeatures, Cstring, ())), \",\"))'"
	    ]
	}
    ]
};

[
    Pipeline("linux", "arm",   "1.3"),
    Pipeline("linux", "arm",   "1.4"),
    Pipeline("linux", "arm",   "1.5"),
    Pipeline("linux", "arm64", "1.3"),
    Pipeline("linux", "arm64", "1.4"),
    Pipeline("linux", "arm64", "1.5")
]
