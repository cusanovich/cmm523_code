#!/bin/bash
# ============================================================================
# CMM 523 -- launch RStudio Server from inside the course container
#
# HOW TO USE:
#   1. In Open OnDemand, start a DESKTOP session (not the RStudio app).
#   2. Open a terminal inside that desktop session.
#   3. Run this script:   bash rserver_cmm523.sh
#   4. It will print a URL. Open Firefox in the desktop session and go there.
#   5. Log in with your NetID and the password printed below.
#   6. Ctrl-C in this terminal when you're done to shut the server down.
# ============================================================================

# --- CHANGE ME -------------------------------------------------------------
PASSWORD="changeme"
CONTAINER=/xdisk/darrenc/cmm_523/containers/cmm523_base.sif

# --- paths -----------------------------------------------------------------
NETID=$(whoami)
WD=/xdisk/darrenc/cmm_523/${NETID}/rstudio
export HUB_ROOT=/xdisk/darrenc/cmm_523/hub_cache

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
    -B /xdisk \
    --env HUB_ROOT="${HUB_ROOT}" \
    "${CONTAINER}" \
    rserver --auth-none=0 \
            --auth-pam-helper-path=pam-helper \
            --server-user="${NETID}" \
            --www-address=127.0.0.1 \
            --www-port="${PORT}"
