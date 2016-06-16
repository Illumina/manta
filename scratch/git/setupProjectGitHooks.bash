#!/usr/bin/env bash

set -u nounset

rel2abs() {
  (cd $1 && pwd -P)
}

thisDir=$(rel2abs $(dirname $0))
baseDir=$(rel2abs $thisDir/../..)
hookDir=$baseDir/.git/hooks

hook1=$hookDir/prepare-commit-msg

if [ -e $hook1 ]; then
    mv $hook1 $hook1.$$.backup
fi

cat << END > $hook1
#!/usr/bin/env bash
#

# this adds the branch name to the beginning of the commit message
BRANCH_NAME=\$(git rev-parse --abbrev-ref HEAD)
FIRST_LINE=\$(head -1 \$1)
if [ -z \$FIRST_LINE ] ; then
    if [ "X\$BRANCH_NAME" != "Xmaster" ] && [ "X\$BRANCH_NAME" != "Xdevelop" ] ; then
        sed -i "1s/^/\$BRANCH_NAME \n/" \$1
    fi
fi
END

chmod +x $hook1

