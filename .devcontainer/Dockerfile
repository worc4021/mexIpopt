FROM ipopt:release AS Rel 
FROM ipopt:debug AS Deb

FROM matlabbuilder

COPY --from=Rel /build/cache /build/cache/Release
COPY --from=Deb /build/cache /build/cache/Debug
COPY pkgconfig.patch /build/cache/pkgconfig.patch
RUN patch -u /build/cache/Release/lib/pkgconfig/ipopt.pc < /build/cache/pkgconfig.patch \
    && patch -u /build/cache/Release/lib/pkgconfig/coinhsl.pc < /build/cache/pkgconfig.patch \
    && patch -u /build/cache/Release/lib/pkgconfig/coinmumps.pc < /build/cache/pkgconfig.patch \
    && patch -u /build/cache/Release/lib/pkgconfig/coinmetis.pc < /build/cache/pkgconfig.patch \
    && patch -u /build/cache/Debug/lib/pkgconfig/ipopt.pc < /build/cache/pkgconfig.patch \
    && patch -u /build/cache/Debug/lib/pkgconfig/coinhsl.pc < /build/cache/pkgconfig.patch \
    && patch -u /build/cache/Debug/lib/pkgconfig/coinmumps.pc < /build/cache/pkgconfig.patch \
    && patch -u /build/cache/Debug/lib/pkgconfig/coinmetis.pc < /build/cache/pkgconfig.patch


ARG USERNAME=vscode

RUN adduser --shell /bin/bash --disabled-password --gecos "" ${USERNAME} \
    && echo "$USERNAME ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/${USERNAME} \
    && chmod 0440 /etc/sudoers.d/${USERNAME}