PS1="touchiyama:\t:@Mac:\w $ "

alias rm='rmtrash'
alias vi='vim'
alias em='emacs'
alias jupyter='jupyter notebook --notebook-dir=~'

RStudio(){
	open /Applications/RStudio.app
}

qtech(){
	ssh -l touchiyama 192.168.0.2
}

qu(){
	ssh -X wakasugi
}

tenki(){
	curl ja.wttr.in/$1
}

when(){
	cal
}

web(){
	ssh bioinfo@bioreg.kyushu-u.ac.jp@hosting4.cc.kyushu-u.ac.jp
}

export PATH="$HOME/.pyenv/shims:$PATH"
export PATH="/usr/local/sbin:$PATH"
export DISPLAY=:0
