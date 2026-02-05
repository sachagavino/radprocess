#!/usr/bin/env bash
set -e

echo "🔧 Setting up radprocess with uv..."

# Ensure uv exists
if ! command -v uv &> /dev/null; then
    echo "❌ uv not found. Install with:"
    echo "   curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
fi

# Create venv if missing
if [ ! -d ".venv" ]; then
    echo "📦 Creating virtual environment..."
    uv venv
fi

# Activate
echo "🐍 Activating venv..."
source .venv/bin/activate

# Install editable
echo "📥 Installing radprocess..."
uv pip install -e .

# Optional: add alias
RADPROCESS_PATH="$(pwd)"

ALIAS_LINE="alias radprocess='uv run --project $RADPROCESS_PATH radprocess'"

if ! grep -q "alias radprocess=" ~/.bashrc 2>/dev/null; then
    echo "$ALIAS_LINE" >> ~/.bashrc
    echo "➕ Added alias to ~/.bashrc"
else
    echo "ℹ️ Alias already exists"
fi

echo
echo "✅ Setup complete!"
echo
echo "You can now:"
echo "  source .venv/bin/activate"
echo "  radprocess"
echo
echo "Or restart shell and just run:"
echo "  radprocess"
