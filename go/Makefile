
GOBINS=lariat
export CGO_LDFLAGS = -L$(shell pwd)/src/gobwa/bwa -L$(shell pwd)/src/jemalloc
export GOPATH=$(shell pwd)

VERSION=$(shell git describe --tags --always --dirty)

$(GOBINS): src/gobwa/bwa/libbwa.a
	go install -ldflags "-X inference.__VERSION__ '$(VERSION)'" $@

src/gobwa/bwa/libbwa.a:
	make -C src/gobwa/bwa libbwa.a

clean:
	rm -Rf bin/ pkg
	$(MAKE) -C src/gobwa/bwa clean

test:
	cd src/test; go test -v
