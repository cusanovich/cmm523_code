#!/bin/bash
# ============================================================================
# CMM 523 -- launch RStudio Server from inside the course container
#
# HOW TO USE:
#   0. Copy this script somewhere you can edit it:
#        cp /groups/darrenc/cmm_523/recipes/rserver_cmm523.sh .
#   1. In Open OnDemand, start a DESKTOP session (not the RStudio app).
#   2. Open a terminal inside that desktop session.
#   3. Run this script:   bash rserver_cmm523.sh
#   4. It will print a URL. Open Firefox in the desktop session and go there.
#   5. Log in with your NetID and the password the script prints.
#   6. Ctrl-C in this terminal when you're done to shut the server down.
# ============================================================================

# --- password ---------------------------------------------------------------
# Generated fresh each launch and printed below. There is nothing to edit and
# nothing to remember: read it off the terminal, type it into the browser.
#
# It is not decoration. RStudio Server listens on this node, and anyone else
# logged into the same compute node could otherwise open your session.
#
# To set your own instead:  RSPASS=mypassword bash rserver_cmm523.sh
PASSWORD=${RSPASS:-$(tr -dc 'A-Za-z0-9' < /dev/urandom | head -c 12)}

# Which container to run RStudio from. Pass YOUR image as the first argument:
#
#   bash rserver_cmm523.sh /xdisk/darrenc/cmm_523/[netid]/week2.sif
#
# With no argument you get the course base image, which deliberately has almost
# nothing installed in it. RStudio will start, but library(Seurat) will fail
# with "there is no package called 'Seurat'". That is the intended behaviour:
# if you forget to name your own container, you should find out immediately
# rather than get plausible-looking results from an environment you did not
# build.
#
# Update the path each week as you add layers -- week2.sif, week3.sif, and so
# on. If your own build is broken and you are stuck, ask; there is a complete
# image available as a fallback.
CONTAINER=${1:-/groups/darrenc/cmm_523/containers/cmm523_base.sif}

# --- paths -----------------------------------------------------------------
NETID=$(whoami)
WD=/xdisk/darrenc/cmm_523/${NETID}/rstudio
export HUB_ROOT=/groups/darrenc/cmm_523/references/hub_cache

# Derive a port from the user id so two people on the same node don't collide
# on 8787. Ports 1024-65535; this keeps us well inside the ephemeral range.
PORT=$(( 8800 + (UID % 1000) ))

if [ ! -f "${CONTAINER}" ]; then
    echo "ERROR: container not found: ${CONTAINER}"
    exit 1
fi

# --- scratch dirs RStudio Server needs to be able to write -----------------
TMPDIR=${WD}/rstudio-tmp
mkdir -p "${TMPDIR}/tmp/rstudio-server" "${TMPDIR}/var/lib" "${TMPDIR}/var/run"

if [ ! -f "${TMPDIR}/tmp/rstudio-server/secure-cookie-key" ]; then
    uuidgen > "${TMPDIR}/tmp/rstudio-server/secure-cookie-key"
    chmod 600 "${TMPDIR}/tmp/rstudio-server/secure-cookie-key"
fi

# --- clear any previous server ---------------------------------------------
# rserver and rsession are separate processes. If an old rsession survives, the
# port stays bound ("Address already in use") AND RStudio resumes the suspended
# session, so config changes appear not to take effect.
# -x matches the process NAME exactly. Do NOT use -f here: -f searches the
# whole command line, which includes "bash rserver_cmm523.sh" -- so the script
# finds itself, kills itself, and prints "Killed" before doing anything.
pkill -9 -u "$USER" -x rserver  2>/dev/null
pkill -9 -u "$USER" -x rsession 2>/dev/null
sleep 2

echo "=================================================="
echo " Open Firefox in this desktop session and go to:"
echo ""
echo "     localhost:${PORT}"
echo ""
echo " Username: ${NETID}"
echo " Password: ${PASSWORD}"
echo ""
echo " Leave this terminal open. Ctrl-C here to stop."
echo "=================================================="

PASSWORD="${PASSWORD}" apptainer exec \
    -B "${TMPDIR}/var/lib:/var/lib/rstudio-server" \
    -B "${TMPDIR}/var/run:/var/run/rstudio-server" \
    -B "${TMPDIR}/tmp:/tmp" \
    -B /xdisk -B /groups \
    --env HUB_ROOT="${HUB_ROOT}" \
    "${CONTAINER}" \
    rserver --auth-none=0 \
            --auth-pam-helper-path=pam-helper \
            --server-user="${NETID}" \
            --www-address=127.0.0.1 \
            --www-port="${PORT}"
