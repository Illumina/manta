
class Manta < Formula
  desc "Structural variant and indel caller for mapped sequencing data"
  homepage "https://github.com/Illumina/manta"
  url "https://github.com/Illumina/manta/releases/download/v0.29.6/manta-0.29.6.release_src.tar.bz2"
  version "0.29.6.release_src"
  sha256 "f06fcb33290f78924c0d87e5a0e3b70150e828246e218fe3869a49684fc9dc16"

  needs :cxx11
  depends_on "boost"
  depends_on "cmake" => :build
  depends_on "zlib"

  def install
    ENV.cxx11
    mkdir "build" do
      system "../configure", "--prefix=#{prefix}"
      system "make", "install"
    end
  end

  test do
    system "python", "${bin}/runMantaWorkflowDemo.py"
  end
end

