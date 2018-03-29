
GOBINS=lariat
export CGO_LDFLAGS = -L$(shell pwd)/src/gobwa/bwa -L$(shell pwd)/src/jemalloc
export GOPATH=$(shell pwd)

VERSION=$(shell git describe --tags --always --dirty)

GO_VERSION=$(strip $(shell go version | sed 's/.*go\([0-9]*\.[0-9]*\).*/\1/'))

# Older versions of Go use the "-X foo bar" syntax.  Newer versions either warn
# or error on that syntax, and use "-X foo=bar" instead.
LINK_SEPERATOR=$(if $(filter 1.5, $(word 1, $(sort 1.5 $(GO_VERSION)))),=, )

$(GOBINS): src/gobwa/bwa/libbwa.a
	go install -ldflags "-X inference.__VERSION__$(LINK_SEPERATOR)'$(VERSION)'" $@

src/gobwa/bwa/libbwa.a:
	make -C src/gobwa/bwa libbwa.a

clean:
	rm -Rf bin/ pkg
	$(MAKE) -C src/gobwa/bwa clean

test:
	cd src/test; go test -v
