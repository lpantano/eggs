require 'formula'

class RnaStar < Formula
  homepage 'https://code.google.com/p/rna-star'
  version '2.4.0j'
  url 'https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.zip'
  sha1 '477238d535c88ad2b8cede69548a21f361d2f39b'

  def install
    if OS.mac?
      bin.install 'bin/MacOSX_x86_64/STAR' => 'STAR'
    else
      bin.install 'bin/Linux_x86_64_static/STAR' => 'STAR'
    end
    bin.install_symlink '../libexec/wrappers/trans-abyss'
  end

end
