#!/bin/bash

echo "🚀 Setting up NextClade environment..."

# Update package lists
sudo apt-get update

# Install NextClade
echo "📦 Installing NextClade..."
curl -fsSL "https://github.com/nextstrain/nextclade/releases/latest/download/nextalign-x86_64-unknown-linux-gnu" -o /tmp/nextalign
curl -fsSL "https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu" -o /tmp/nextclade

# Make executables
sudo chmod +x /tmp/nextalign /tmp/nextclade
sudo mv /tmp/nextalign /usr/local/bin/nextalign
sudo mv /tmp/nextclade /usr/local/bin/nextclade

# Verify installation
echo "✅ Verifying NextClade installation..."
nextclade --version
nextalign --version


echo "✨ Setup complete! NextClade is ready to use."
echo "Try running: nextclade --help"
